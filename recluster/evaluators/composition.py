"""composition_shift: does the population's share differ between conditions, across replicates?

The unit of replication is the sample (donor x condition), never the cell. A naive cell-level
chi-square is computed only to show how misleading it is.
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
) -> pd.DataFrame:
    """One row per sample: n cells in run, k cells in population, proportion, condition."""
    s_run = sample.loc[run_cells].astype(str)
    n = s_run.value_counts()
    k = sample.loc[cells_in_population].astype(str).value_counts().reindex(n.index, fill_value=0)
    cond = condition.loc[run_cells].astype(str).groupby(s_run.to_numpy()).first().reindex(n.index)
    df = pd.DataFrame({"n": n, "k": k, "condition": cond})
    df["proportion"] = df["k"] / df["n"]
    df["logit"] = np.log((df["k"] + 0.5) / (df["n"] - df["k"] + 0.5))
    return df.sort_index()


def evaluate_composition(
    cells_in_population: pd.Index,
    run_cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    comparison: Comparison,
    thr: Thresholds,
) -> dict[str, object]:
    table = per_sample_proportions(cells_in_population, run_cells, sample, condition)
    a = table[table["condition"] == comparison.group_a]
    b = table[table["condition"] == comparison.group_b]
    warnings: list[str] = []

    effect_logodds = float(a["logit"].mean() - b["logit"].mean()) if len(a) and len(b) else math.nan
    mean_a, mean_b = (
        float(a["proportion"].mean()) if len(a) else math.nan,
        (float(b["proportion"].mean()) if len(b) else math.nan),
    )
    ratio = mean_a / mean_b if mean_b and not math.isnan(mean_b) and mean_b > 0 else math.nan

    if len(a) >= 2 and len(b) >= 2:
        t_p = float(stats.ttest_ind(a["logit"], b["logit"], equal_var=False).pvalue)
        mw_p = float(
            stats.mannwhitneyu(a["proportion"], b["proportion"], alternative="two-sided").pvalue
        )
    else:
        t_p = mw_p = math.nan
        warnings.append(
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
    holds = bool(not math.isnan(decisive_p) and decisive_p < thr.alpha and direction_ok)
    holds_cell = bool(not math.isnan(chi_p) and chi_p < thr.alpha and direction_ok)

    return {
        "holds": holds,
        "holds_cell_level": holds_cell,
        "replicate_ok": holds,
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
