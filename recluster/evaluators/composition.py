"""composition_shift: does the population's share differ between conditions, across replicates?

The unit of replication is the sample (donor x condition), never the cell. When the same
donors appear in both groups the test is paired on the donor. A run can only test the claim
when the population is present in enough samples of both groups: when the matching finds it in
one condition only, the tool cannot tell a biological absence from the clustering having lost
that condition's version of the population, and the first benchmark showed exactly that
artifact producing a huge spurious shift. Such runs are reported as not testable. A naive
cell-level chi-square is computed only to show how misleading it is.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

from recluster.config import Thresholds
from recluster.spec import Comparison


def per_sample_proportions(
    cells_in_population: pd.Index,
    run_cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    pair: pd.Series | None = None,
) -> pd.DataFrame:
    """One row per sample: n cells in run, k cells in population, proportion, condition, pair."""
    s_run = sample.loc[run_cells].astype(str)
    n = s_run.value_counts()
    k = sample.loc[cells_in_population].astype(str).value_counts().reindex(n.index, fill_value=0)
    by_sample = s_run.to_numpy()
    cond = condition.loc[run_cells].astype(str).groupby(by_sample).first().reindex(n.index)
    df = pd.DataFrame({"n": n, "k": k, "condition": cond})
    if pair is not None:
        # the pair key covers the two compared groups only; other cells get no pair id
        by_cell = pair.astype(str).reindex(run_cells)
        df["pair"] = by_cell.groupby(by_sample).first().reindex(n.index)
    df["proportion"] = df["k"] / df["n"]
    df["logit"] = np.log((df["k"] + 0.5) / (df["n"] - df["k"] + 0.5))
    return df.sort_index()


def group_presence(
    table: pd.DataFrame, comparison: Comparison, thr: Thresholds
) -> tuple[bool, list[str]]:
    """Is the population present (>= min cells) in enough samples of both groups?

    Absence from one group is not accepted as evidence of a shift: in a re-clustering the
    population can vanish from one condition because the clusters of that condition merged
    elsewhere, and a per-sample test would then report a large, spurious difference (the first
    benchmark's only composition false alarm arose this way). The run is reported as not
    testable instead; the verdict counts such runs as failures and says how many there were.
    """
    reasons = []
    for grp in (comparison.group_a, comparison.group_b):
        sub = table[table["condition"] == grp]
        present = int((sub["k"] >= thr.min_cells_per_sample).sum())
        if present < thr.min_samples_present:
            reasons.append(
                f"population present in only {present} sample(s) of group {grp!r} "
                f"(need {thr.min_samples_present}); the run cannot test this claim"
            )
    return (not reasons), reasons


def evaluate_composition(
    cells_in_population: pd.Index,
    run_cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    comparison: Comparison,
    thr: Thresholds,
    *,
    pair: pd.Series | None = None,
) -> dict[str, object]:
    table = per_sample_proportions(cells_in_population, run_cells, sample, condition, pair)
    a = table[table["condition"] == comparison.group_a]
    b = table[table["condition"] == comparison.group_b]
    warnings: list[str] = []
    _, reasons = group_presence(table, comparison, thr)

    effect_logodds = float(a["logit"].mean() - b["logit"].mean()) if len(a) and len(b) else math.nan
    mean_a = float(a["proportion"].mean()) if len(a) else math.nan
    mean_b = float(b["proportion"].mean()) if len(b) else math.nan
    ratio = mean_a / mean_b if mean_b and not math.isnan(mean_b) and mean_b > 0 else math.nan

    paired = pair is not None and "pair" in table.columns
    t_p = mw_p = math.nan
    design = "unpaired"
    design_note: str | None = None
    if paired:
        diff = _paired_differences(a, b)
        if len(diff) >= 2:
            design = "paired"
            t_p = float(stats.ttest_rel(diff["logit_a"], diff["logit_b"]).pvalue)
            d = diff["logit_a"] - diff["logit_b"]
            mw_p = _signed_rank_p(d.to_numpy())
            effect_logodds = float(d.mean())
        else:
            design_note = "fewer than 2 complete pairs; falling back to the unpaired test"
            warnings.append(design_note)
            paired = False
    if not paired:
        if len(a) >= 2 and len(b) >= 2:
            t_p = float(stats.ttest_ind(a["logit"], b["logit"], equal_var=False).pvalue)
            mw_p = float(
                stats.mannwhitneyu(a["proportion"], b["proportion"], alternative="two-sided").pvalue
            )
        else:
            reasons.append(
                f"fewer than 2 replicates per group ({len(a)} vs {len(b)}); "
                "replicate-aware test undefined"
            )

    # naive cell-level contrast: chi-square on pooled cells (pseudoreplication)
    ka, na = int(a["k"].sum()), int(a["n"].sum())
    kb, nb = int(b["k"].sum()), int(b["n"].sum())
    if na and nb and (ka + kb) > 0 and (na - ka + nb - kb) > 0:
        chi_p = float(stats.chi2_contingency([[ka, na - ka], [kb, nb - kb]], correction=True)[1])
    else:
        chi_p = math.nan

    direction_ok = _direction_ok(effect_logodds, comparison.direction)
    decisive_p = t_p if thr.composition_test == "logit_t" else mw_p
    evaluable = not reasons
    holds = bool(
        evaluable and not math.isnan(decisive_p) and decisive_p < thr.alpha and direction_ok
    )
    holds_cell = bool(not math.isnan(chi_p) and chi_p < thr.alpha and direction_ok)
    warnings.extend(reasons)

    return {
        "holds": holds,
        "holds_cell_level": holds_cell,
        "replicate_ok": holds,
        "evaluable": evaluable,
        "reasons": reasons,
        "design": design,
        "design_note": design_note,
        "direction_ok": direction_ok,
        "effect_log_odds_a_minus_b": _r(effect_logodds),
        "mean_proportion_a": _r(mean_a),
        "mean_proportion_b": _r(mean_b),
        "proportion_ratio_a_over_b": _r(ratio),
        "p_logit_ttest": _r(t_p),
        "p_mannwhitney": _r(mw_p),
        "p_cell_level_chi2": _r(chi_p),
        "n_samples_a": int(len(a)),
        "n_samples_b": int(len(b)),
        "n_cells_a": na,
        "n_cells_b": nb,
        "per_sample": table.reset_index()
        .rename(columns={"index": "sample"})
        .to_dict(orient="records"),
        "warnings": warnings,
        "key_stats": {"effect_log_odds": _r(effect_logodds), "p_replicate": _r(decisive_p)},
    }


def _paired_differences(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """One row per pair level present in both groups (first sample per level on each side)."""
    pa = a.groupby("pair")["logit"].first().rename("logit_a")
    pb = b.groupby("pair")["logit"].first().rename("logit_b")
    return pd.concat([pa, pb], axis=1, join="inner")


def _signed_rank_p(d: np.ndarray) -> float:
    d = d[d != 0]
    if len(d) < 2:
        return math.nan
    try:
        return float(stats.wilcoxon(d, alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def _direction_ok(effect: float, direction: str) -> bool:
    if math.isnan(effect):
        return False
    if direction == "higher_in_a":
        return effect > 0
    if direction == "lower_in_a":
        return effect < 0
    return True


def _r(x: float) -> float | None:
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else float(x)
