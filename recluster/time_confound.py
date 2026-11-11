"""Time-of-day confound check.

Human gene expression changes with time of collection: core clock genes differ between day
and night across many tissues (Wucher et al., PLoS Biol 2023;21(2):e3001986, GTEx, 932
donors). If one condition was sampled mostly in the morning and the other in the afternoon, a
"condition effect" can be a time effect. Nothing in a standard single-cell workflow checks this.

Two modes, both at the level of the biological replicate (never the cell):

Mode A, collection time known. (1) Balance: is time associated with condition across samples?
Each sample's clock time becomes an angle; the groups' mean (cos, sin) vectors are compared with
a permutation test (exact when at most 20,000 label assignments exist, otherwise 5,000 seeded
draws), and the circular mean difference is reported in hours. The groups count as imbalanced
when that test is significant or when the circular means differ by at least
``Thresholds.time_imbalance_hours``: with three or fewer samples per group no permutation test
can reach alpha, and the hours rule keeps a large practical difference from passing as
"balanced". (2) Adjustment: the claim's replicate-level model (the same design the verdict
used) is refitted with cos(2*pi*t/24) and sin(2*pi*t/24) as covariates. The report says whether
the effect survives, how much its estimate shrinks, and whether its direction flips. A design in
which the time covariates cannot be separated from the condition (every control sample at one
hour and every treated sample at another) is reported as such, not as an effect that vanished.

Mode B, collection time unknown. Per-sample expression of a panel of 13 canonical core clock
genes (the transcriptional-translational feedback loop reviewed in Takahashi, Nat Rev Genet
2017;18:164-179; the same genes Wucher et al. 2023 report as day-night variable) is tested
between groups, gene by gene and as a phase-structured index (E-box arm minus activator arm).
This is a red-flag screen: it says whether clock-gene state differs between the groups in the
pattern a collection-time difference produces. It does NOT estimate circadian phase, which is
not reliable from single-cell data without a validated method.

Units. The unit of both permutation tests is the sample, or the donor when samples are nested
in donors that share one condition (and one collection time): two samples of one donor are not
two independent observations of when that donor was sampled.

Rhythmic-gene overlap. When a public list of day-night variable genes is available (the
Wucher et al. 2023 supplementary tables, CC BY 4.0), claimed genes on that list are flagged.
Without a usable list the step is skipped and the report says so.
"""

from __future__ import annotations

import contextlib
import io
import math
import warnings
from dataclasses import asdict, dataclass, field
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from recluster.config import Thresholds
from recluster.evaluators.pseudobulk import pseudobulk_counts
from recluster.hours import parse_hours
from recluster.spec import Comparison

# Core clock panel. Arms follow the canonical loop: the activator arm (BMAL1/CLOCK/NPAS2)
# peaks roughly in antiphase to the E-box-driven repressor/output arm (PER, CRY, REV-ERB, DBP,
# TEF, CIART) (Takahashi 2017, Fig. 1).
CLOCK_PANEL: dict[str, str] = {
    "ARNTL": "activator",  # BMAL1
    "CLOCK": "activator",
    "NPAS2": "activator",
    "PER1": "ebox",
    "PER2": "ebox",
    "PER3": "ebox",
    "CRY1": "ebox",
    "CRY2": "ebox",
    "NR1D1": "ebox",  # REV-ERB alpha
    "NR1D2": "ebox",  # REV-ERB beta
    "DBP": "ebox",
    "TEF": "ebox",
    "CIART": "ebox",
}
# CRY1 is tested gene by gene but left out of the clock-state index: its transcript peaks
# several hours after the other E-box targets, between the two arms (Ukai-Tadenuma et al.,
# Cell 2011;144:268-281), so it would blur an index built on antiphase.
INDEX_EXCLUDED: frozenset[str] = frozenset({"CRY1"})
PANEL_CITATION = (
    "Takahashi JS, Nat Rev Genet 2017;18:164-179; Wucher V et al., PLoS Biol 2023;21(2):e3001986"
)

NO_CONFOUND = "NO TIME CONFOUND DETECTED"
POSSIBLE = "POSSIBLE TIME CONFOUND"
NOT_ASSESSABLE = "NOT ASSESSABLE"
WARNING_TEXT = (
    "clock-gene state differs between groups, so a time-of-collection or circadian confound "
    "is possible"
)
MAX_EXACT_PERMUTATIONS = 20000
MONTE_CARLO_DRAWS = 5000
MAX_TIME_SPREAD_HOURS = 0.5  # cells of one sample (or samples of one donor) must agree this well


# --------------------------------------------------------------------------- parsing


def _circular_mean_hours(theta: np.ndarray) -> float:
    v = _mean_vector(theta)
    hours = float((math.atan2(v[1], v[0]) % (2 * math.pi)) / (2 * math.pi) * 24)
    return 0.0 if math.isclose(hours, 24.0, abs_tol=1e-9) else hours % 24.0


def _circular_mean_of_hours(hours: pd.Series | np.ndarray) -> float:
    v = np.asarray(pd.Series(hours).dropna(), dtype=float)
    if v.size == 0:
        return math.nan
    return _circular_mean_hours(v / 24.0 * 2 * math.pi)


def _circular_range_hours(hours: pd.Series | np.ndarray) -> float:
    """Largest clock distance between any two of the hours (so 23:30 and 00:30 are 1 h apart)."""
    v = np.asarray(pd.Series(hours).dropna(), dtype=float)
    if v.size < 2:
        return 0.0
    d = np.abs(v[:, None] - v[None, :]) % 24.0
    return float(np.minimum(d, 24.0 - d).max())


def sample_times(time_per_cell: pd.Series, sample: pd.Series) -> tuple[pd.Series, list[str]]:
    """One collection time per sample: the circular mean of its cells' times.

    The per-cell values are expected to agree within a sample; a spread above half an hour is
    reported. The mean is circular so that a sample straddling midnight is not placed at noon.
    """
    hours = parse_hours(time_per_cell)
    df = pd.DataFrame({"t": hours.to_numpy(), "s": sample.loc[hours.index].astype(str).to_numpy()})
    warns = []
    per = df.groupby("s")["t"].agg(
        mean=_circular_mean_of_hours,
        spread=_circular_range_hours,
        na=lambda x: float(x.isna().mean()),
    )
    bad = per[per["spread"] > MAX_TIME_SPREAD_HOURS]
    if len(bad):
        warns.append(
            f"{len(bad)} sample(s) carry more than one collection time; the circular mean was used"
        )
    if (per["na"] > 0).any():
        warns.append(f"{int((per['na'] > 0).sum())} sample(s) have unparseable collection times")
    return per["mean"].astype(float), warns


# ----------------------------------------------------------------------------- units


def exchangeable_units(
    obs: pd.DataFrame,
    cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    group_a: str,
    group_b: str,
    sample_columns: list[str] | None,
    condition_column: str | None,
    hours_by_sample: pd.Series | None = None,
) -> tuple[pd.Series, str | None, list[str]]:
    """The units of the sample-level permutation tests: one id per cell.

    The unit is the sample unless a sample-key column other than the condition column groups
    several samples per level while the condition (and the collection time, when given) is
    constant within every level and each group holds at least two levels; the coarsest such
    column wins. Example: ``sample = [donor, label]`` with a condition assigned per donor makes
    the donor the unit, because a donor's two samples are not two independent observations of
    its condition or of when it was sampled. Returns the unit id per cell (cells of the two
    compared groups only), the nesting column or None, and notes for the report.
    """
    cond = condition.loc[cells].astype(str)
    cells = cells[cond.isin([group_a, group_b]).to_numpy()]
    cond = cond.loc[cells]
    samp = sample.loc[cells].astype(str)
    best: tuple[int, str, pd.Series] | None = None
    for col in sample_columns or []:
        if col == condition_column or col not in obs.columns:
            continue
        level = obs.loc[cells, col].astype(str)
        frame = pd.DataFrame(
            {"level": level.to_numpy(), "cond": cond.to_numpy(), "sample": samp.to_numpy()}
        )
        by_level = frame.groupby("level")
        if (by_level["cond"].nunique() > 1).any():
            continue  # condition varies within a level: the level is not the unit
        if int(by_level["sample"].nunique().max()) < 2:
            continue  # no coarser than the sample itself
        cond_of_level = by_level["cond"].first()
        if (cond_of_level == group_a).sum() < 2 or (cond_of_level == group_b).sum() < 2:
            continue
        if hours_by_sample is not None:
            per_sample = frame.drop_duplicates("sample")
            hrs = hours_by_sample.reindex(per_sample["sample"].to_numpy()).to_numpy()
            spread = (
                pd.DataFrame({"level": per_sample["level"].to_numpy(), "h": hrs})
                .groupby("level")["h"]
                .agg(_circular_range_hours)
            )
            if (spread > MAX_TIME_SPREAD_HOURS).any():
                continue  # time varies within a level
        n_levels = int(len(cond_of_level))
        if best is None or n_levels < best[0]:
            best = (n_levels, col, level.rename("unit"))
    if best is None:
        return samp.rename("unit"), None, []
    n_levels, col, unit = best
    note = (
        f"samples are nested in `{col}` ({int(samp.nunique())} samples in {n_levels} levels, "
        "each level in one condition"
        + (" with one collection time" if hours_by_sample is not None else "")
        + f"): the {n_levels} levels are the units of the sample-level tests"
    )
    return unit, col, [note]


def _unit_of_sample(unit: pd.Series, sample: pd.Series) -> pd.Series:
    """Map each sample to its unit (units never split a sample)."""
    s = sample.loc[unit.index].astype(str).to_numpy()
    return unit.groupby(s).first()


# --------------------------------------------------------------------------- balance


def _mean_vector(theta: np.ndarray) -> np.ndarray:
    return np.array([np.cos(theta).mean(), np.sin(theta).mean()])


def _permutation_p(
    values: np.ndarray, is_a: np.ndarray, statistic: Any, *, seed: int = 0
) -> dict[str, Any]:
    """Permutation p for a label-symmetric statistic of (values, is_a).

    Exact over every label assignment when there are at most ``MAX_EXACT_PERMUTATIONS``,
    otherwise ``MONTE_CARLO_DRAWS`` seeded draws. ``min_attainable_p`` is the smallest p the
    test can produce (two assignments tie at the maximum when the groups are equal in size, as
    the statistic does not change under swapping the labels); a test whose minimum is at or
    above alpha cannot reject, whatever the data.
    """
    n, k = len(is_a), int(is_a.sum())
    observed = statistic(values, is_a)
    n_total = math.comb(n, k)
    count = 0
    if n_total <= MAX_EXACT_PERMUTATIONS:
        for idx in combinations(range(n), k):
            mask = np.zeros(n, dtype=bool)
            mask[list(idx)] = True
            if statistic(values, mask) >= observed - 1e-12:
                count += 1
        return {
            "p": count / n_total,
            "n_permutations": n_total,
            "exact": True,
            "min_attainable_p": (2 if 2 * k == n else 1) / n_total,
        }
    rng = np.random.default_rng(seed)
    for _ in range(MONTE_CARLO_DRAWS):
        mask = np.zeros(n, dtype=bool)
        mask[rng.choice(n, size=k, replace=False)] = True
        if statistic(values, mask) >= observed - 1e-12:
            count += 1
    return {
        "p": (count + 1) / (MONTE_CARLO_DRAWS + 1),
        "n_permutations": MONTE_CARLO_DRAWS,
        "exact": False,
        "min_attainable_p": 1 / (MONTE_CARLO_DRAWS + 1),
    }


def balance_check(
    hours_by_unit: pd.Series,
    group_by_unit: pd.Series,
    group_a: str,
    group_b: str,
    thr: Thresholds | None = None,
    *,
    unit_column: str | None = None,
) -> dict[str, Any]:
    """Is collection time associated with the condition across the units?

    Imbalanced when the permutation test on the groups' mean (cos, sin) vectors is significant
    or when the circular mean times differ by at least ``thr.time_imbalance_hours``. Both are
    reported, with the smallest p the permutation test could have produced.
    """
    thr = thr or Thresholds()
    hours = hours_by_unit.astype(float)
    groups = group_by_unit.reindex(hours.index)
    keep = groups.isin([group_a, group_b]) & hours.notna()
    h, g = hours[keep], groups[keep]
    theta = (h.to_numpy() / 24.0) * 2 * math.pi
    is_a = (g == group_a).to_numpy()
    n_a, n_b = int(is_a.sum()), int((~is_a).sum())
    out: dict[str, Any] = {"unit_column": unit_column, "n_units_a": n_a, "n_units_b": n_b}
    if n_a < 2 or n_b < 2:
        out.update(
            {
                "assessable": False,
                "reason": f"fewer than 2 units with a collection time per group ({n_a} vs {n_b})",
            }
        )
        return out

    def stat(th: np.ndarray, mask: np.ndarray) -> float:
        return float(np.linalg.norm(_mean_vector(th[mask]) - _mean_vector(th[~mask])))

    perm = _permutation_p(theta, is_a, stat)
    mean_a, mean_b = _circular_mean_hours(theta[is_a]), _circular_mean_hours(theta[~is_a])
    diff = float(((mean_a - mean_b + 12) % 24) - 12)
    by_test = bool(perm["p"] < thr.alpha)
    by_hours = bool(abs(diff) >= thr.time_imbalance_hours)
    if by_test and by_hours:
        basis = f"permutation test and mean difference of at least {thr.time_imbalance_hours:g} h"
    elif by_test:
        basis = "permutation test"
    elif by_hours:
        basis = f"mean difference of at least {thr.time_imbalance_hours:g} h"
    else:
        basis = "neither"
    out.update(
        {
            "assessable": True,
            "circular_mean_hours_a": round(mean_a, 2),
            "circular_mean_hours_b": round(mean_b, 2),
            "circular_mean_difference_hours": round(diff, 2),
            "mean_vector_distance": round(stat(theta, is_a), 4),
            "p_permutation": perm["p"],
            "n_permutations": perm["n_permutations"],
            "exact": perm["exact"],
            "min_attainable_p": perm["min_attainable_p"],
            "test_can_reject": bool(perm["min_attainable_p"] < thr.alpha),
            "imbalanced": by_test or by_hours,
            "imbalance_basis": basis,
            "p_cos_ttest": _welch_p(np.cos(theta[is_a]), np.cos(theta[~is_a])),
            "p_sin_ttest": _welch_p(np.sin(theta[is_a]), np.sin(theta[~is_a])),
            "times_by_unit": {str(u): round(float(t), 2) for u, t in h.items()},
        }
    )
    return out


def _welch_p(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2 or (a.std() == 0 and b.std() == 0):
        return math.nan
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


# ------------------------------------------------------------------------ adjustment


def _pair_levels(
    pair: pd.Series | None, cells: pd.Index, by_sample: np.ndarray, index: pd.Index
) -> pd.Series | None:
    """One pair label per pseudobulk sample, or None when a sample has none."""
    if pair is None:
        return None
    per_cell = pair.reindex(cells)
    if per_cell.isna().any():
        return None
    per_sample = per_cell.astype(str).groupby(by_sample).first().reindex(index)
    if per_sample.isna().any():
        return None
    levels = {p: f"p{i}" for i, p in enumerate(sorted(per_sample.unique()))}
    return per_sample.map(levels)


def _design_columns(
    cond: np.ndarray, theta: np.ndarray, pairs: pd.Series | None
) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    """Numeric design matrices (base, adjusted) with the condition column last.

    Returns the blocking terms ("pair" when pair dummies are present) and the time terms kept:
    a time column that adds no rank to the intercept, the pair dummies and the earlier time
    column is dropped (with two distinct clock times sin is a linear function of cos, and with
    one time per pair both are absorbed by the pairing), so that redundancy in the time basis is
    not mistaken for confounding of the condition.
    """
    cols: list[np.ndarray] = [np.ones_like(cond)]
    terms: list[str] = []
    if pairs is not None:
        dummies = pd.get_dummies(pairs, drop_first=True).to_numpy(dtype=float)
        if dummies.shape[1]:
            cols.append(dummies)
            terms.append("pair")
    block = np.column_stack(cols)
    time_cols = {"cos_t": np.cos(theta), "sin_t": np.sin(theta)}
    kept: list[str] = []
    for name, col in time_cols.items():
        trial = np.column_stack([block] + [time_cols[k] for k in kept] + [col])
        if np.linalg.matrix_rank(trial) == trial.shape[1]:
            kept.append(name)
    base = np.column_stack([block, cond])
    adjusted = np.column_stack([block] + [time_cols[k] for k in kept] + [cond])
    return base, adjusted, terms, kept


def _not_estimable_reason(
    x_base: np.ndarray, x_adj: np.ndarray, time_terms: list[str]
) -> str | None:
    """Why the condition coefficient of the adjusted model cannot be estimated, or None."""
    n = x_adj.shape[0]
    if n <= x_adj.shape[1] or n <= x_base.shape[1]:
        return (
            f"{n} replicates with a collection time are too few for the adjusted design "
            f"({x_adj.shape[1]} terms)"
        )
    if not time_terms:
        return "collection time adds nothing to the model (one distinct time per pair or overall)"
    for x in (x_base, x_adj):
        if np.linalg.matrix_rank(x) < x.shape[1]:
            return (
                "collection time is collinear with the condition (or with the pairing): the "
                "effect cannot be separated from collection time"
            )
    return None


def adjusted_pseudobulk_de(
    counts: AnnData,
    cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    hours_by_sample: pd.Series,
    comparison: Comparison,
    genes: list[str],
    thr: Thresholds,
    *,
    n_cpus: int,
    pair: pd.Series | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Pseudobulk DESeq2 for the claimed genes with and without cos/sin time covariates.

    The base design is the one the verdict used (``~condition``, or ``~pair + condition`` when
    every pair has a pseudobulk in both groups); the adjusted design adds ``cos_t`` and
    ``sin_t``. The refit is skipped, with the reason recorded in the returned info, when there
    are too few pseudobulks for the adjusted design or its columns are collinear (every control
    sample at one hour and every treated sample at another, or one time per pair): the effect
    then cannot be separated from time, which is a finding in itself. Outlier (Cook's) filtering
    is off so that every claimed gene keeps a p-value in both fits.
    """
    pb, _ = pseudobulk_counts(counts, cells, sample, min_cells=thr.min_cells_per_sample)
    by_sample = sample.loc[cells].astype(str).to_numpy()
    cond = condition.loc[cells].astype(str).groupby(by_sample).first().reindex(pb.index)
    hours = hours_by_sample.reindex(pb.index).astype(float)
    keep = cond.isin([comparison.group_a, comparison.group_b]) & hours.notna()
    pb, cond, hours = pb.loc[keep], cond.loc[keep], hours.loc[keep]
    n_a = int((cond == comparison.group_a).sum())
    n_b = int((cond == comparison.group_b).sum())
    empty = pd.DataFrame(
        {
            "lfc_unadjusted": math.nan,
            "p_unadjusted": math.nan,
            "lfc_adjusted": math.nan,
            "p_adjusted": math.nan,
        },
        index=pd.Index(genes, name="gene"),
    )
    info: dict[str, Any] = {
        "design": "~condition",
        "estimable": False,
        "reason": None,
        "n_samples_a": n_a,
        "n_samples_b": n_b,
    }
    if n_a < thr.min_samples_present or n_b < thr.min_samples_present:
        info["reason"] = (
            f"population present in {n_a} vs {n_b} pseudobulk replicates with a collection time "
            f"(need {thr.min_samples_present} per group)"
        )
        return empty, info
    totals = pb.sum(axis=0)
    pb = pb.loc[:, (totals >= thr.pseudobulk_min_total_counts) | pb.columns.isin(genes)]
    is_a = (cond == comparison.group_a).to_numpy()
    theta = 2 * math.pi * hours.to_numpy() / 24.0
    meta = pd.DataFrame(
        {
            "condition": np.where(is_a, "grpA", "grpB"),
            "cos_t": np.cos(theta),
            "sin_t": np.sin(theta),
        },
        index=pb.index,
    )
    pairs = _pair_levels(pair, cells, by_sample, pb.index)
    if pairs is not None:
        complete = pd.Series(meta["condition"].to_numpy()).groupby(pairs.to_numpy()).nunique()
        if not ((complete == 2).all() and len(complete) >= 2):
            pairs = None
    if pairs is not None:
        meta["pair"] = pairs.to_numpy()
    x_base, x_adj, terms, time_terms = _design_columns(is_a.astype(float), theta, pairs)
    info["design"] = "~" + " + ".join(terms + ["condition"])
    info["time_terms"] = time_terms
    reason = _not_estimable_reason(x_base, x_adj, time_terms)
    if reason is not None:
        info["reason"] = reason
        return empty, info
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    def fit(design: str) -> pd.DataFrame:
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            dds = DeseqDataSet(counts=pb, metadata=meta, design=design, quiet=True, n_cpus=n_cpus)
            dds.deseq2()
            st = DeseqStats(
                dds,
                contrast=["condition", "grpA", "grpB"],
                cooks_filter=False,
                independent_filter=False,
                quiet=True,
                n_cpus=n_cpus,
            )
            st.summary()
            return st.results_df.reindex(genes)

    try:
        base = fit(info["design"])
        adj = fit("~" + " + ".join(terms + time_terms + ["condition"]))
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        info["reason"] = f"the refit failed ({exc.__class__.__name__}: {exc})"
        return empty, info
    info["estimable"] = True
    out = pd.DataFrame(
        {
            "lfc_unadjusted": base["log2FoldChange"],
            "p_unadjusted": base["pvalue"],
            "lfc_adjusted": adj["log2FoldChange"],
            "p_adjusted": adj["pvalue"],
        }
    )
    out.index.name = "gene"
    return out, info


def adjusted_composition(
    table: pd.DataFrame, hours_by_sample: pd.Series, group_a: str, group_b: str
) -> dict[str, Any]:
    """OLS of the logit proportion on condition, with and without cos/sin of collection time.

    Samples of other conditions are dropped first. When the table carries a complete pairing
    (every pair in both groups) the pair enters as fixed effects, as in the verdict's paired
    test. The condition coefficient and its p are returned for both fits; a rank-deficient
    design (times collinear with the condition or constant within pairs) is reported as not
    estimable instead of as a vanished effect.
    """
    sub = table[table["condition"].isin([group_a, group_b])]
    t = hours_by_sample.reindex(sub.index).astype(float)
    sub = sub.loc[t.notna().to_numpy()]
    t = t.loc[sub.index]
    y = sub["logit"].to_numpy(dtype=float)
    cond = (sub["condition"] == group_a).to_numpy(dtype=float)
    theta = 2 * math.pi * t.to_numpy() / 24.0
    pairs: pd.Series | None = None
    if "pair" in sub.columns and sub["pair"].notna().all() and len(sub):
        complete = sub.groupby("pair")["condition"].nunique()
        if (complete == 2).all() and len(complete) >= 2:
            pairs = sub["pair"].astype(str)
    x_base, x_adj, terms, time_terms = _design_columns(cond, theta, pairs)
    out: dict[str, Any] = {
        "design": "paired" if "pair" in terms else "unpaired",
        "time_terms": time_terms,
        "n_samples_a": int(cond.sum()),
        "n_samples_b": int(len(cond) - cond.sum()),
        "estimable": False,
        "reason": None,
    }
    reason = _not_estimable_reason(x_base, x_adj, time_terms)
    if reason is not None:
        out["reason"] = reason
        out.update(
            {
                "effect_unadjusted": math.nan,
                "p_unadjusted": math.nan,
                "effect_adjusted": math.nan,
                "p_adjusted": math.nan,
            }
        )
        return out

    def ols(x: np.ndarray) -> tuple[float, float]:
        n, p = x.shape
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        resid = y - x @ beta
        s2 = float(resid @ resid) / (n - p)
        cov = s2 * np.linalg.inv(x.T @ x)
        se = math.sqrt(max(cov[-1, -1], 0.0))
        if se == 0:
            return float(beta[-1]), math.nan
        return float(beta[-1]), float(2 * stats.t.sf(abs(beta[-1] / se), n - p))

    b0, p0 = ols(x_base)
    b1, p1 = ols(x_adj)
    out.update(
        {
            "estimable": True,
            "effect_unadjusted": b0,
            "p_unadjusted": p0,
            "effect_adjusted": b1,
            "p_adjusted": p1,
        }
    )
    return out


# ----------------------------------------------------------------------- clock panel


def panel_expression(
    counts: AnnData, sample: pd.Series, *, min_total_counts: int
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Per-sample log2 CPM of the panel genes from all cells; which genes are usable."""
    present = [g for g in CLOCK_PANEL if g in counts.var_names]
    pb, _ = pseudobulk_counts(counts, counts.obs_names, sample, min_cells=1)
    lib = pb.sum(axis=1)
    sub = pb[present]
    usable = [g for g in present if int(sub[g].sum()) >= min_total_counts]
    cpm = sub.div(lib, axis=0) * 1e6
    return np.log2(cpm + 1.0), sub, usable


def index_genes(usable: list[str]) -> tuple[list[str], list[str]]:
    """Panel genes that enter the clock-state index, by arm."""
    ebox = [g for g in usable if CLOCK_PANEL[g] == "ebox" and g not in INDEX_EXCLUDED]
    act = [g for g in usable if CLOCK_PANEL[g] == "activator" and g not in INDEX_EXCLUDED]
    return ebox, act


def clock_state_index(log_cpm: pd.DataFrame, usable: list[str]) -> pd.Series | None:
    """E-box arm minus activator arm, each the mean of per-gene z-scores across samples."""
    ebox, act = index_genes(usable)
    if not ebox or not act:
        return None
    used = ebox + act
    z = (log_cpm[used] - log_cpm[used].mean()) / log_cpm[used].std(ddof=0).replace(0, 1)
    return (z[ebox].mean(axis=1) - z[act].mean(axis=1)).rename("clock_state_index")


def panel_check(
    counts: AnnData,
    sample: pd.Series,
    group_by_sample: pd.Series,
    group_a: str,
    group_b: str,
    thr: Thresholds,
    *,
    unit_column: str | None = None,
) -> dict[str, Any]:
    """Mode B: does clock-gene state differ between the groups at the sample level?

    ``sample`` is the unit of the test (the sample, or the donor when samples are nested in
    donors of one condition); ``group_by_sample`` gives each unit's condition.
    """
    log_cpm, raw, usable = panel_expression(
        counts, sample, min_total_counts=thr.time_min_gene_counts
    )
    ebox, act = index_genes(usable)
    out: dict[str, Any] = {
        "panel": list(CLOCK_PANEL),
        "panel_citation": PANEL_CITATION,
        "genes_present": list(log_cpm.columns),
        "genes_usable": usable,
        "index_genes": ebox + act,
        "index_excluded": sorted(INDEX_EXCLUDED),
        "min_total_counts": thr.time_min_gene_counts,
        "unit_column": unit_column,
    }
    groups = group_by_sample.reindex(log_cpm.index)
    keep = groups.isin([group_a, group_b])
    log_cpm, groups = log_cpm.loc[keep], groups.loc[keep]
    is_a = (groups == group_a).to_numpy()
    n_a, n_b = int(is_a.sum()), int((~is_a).sum())
    out.update({"n_samples_a": n_a, "n_samples_b": n_b})
    index = clock_state_index(log_cpm, usable)
    if len(usable) < thr.time_min_panel_genes or index is None or n_a < 2 or n_b < 2:
        out["assessable"] = False
        out["reason"] = (
            f"{len(usable)} of {len(CLOCK_PANEL)} panel genes have at least "
            f"{thr.time_min_gene_counts} counts across samples (need {thr.time_min_panel_genes}, "
            "with both clock arms represented) and at least 2 units per group"
        )
        return out

    # per-gene differences across units
    rows = []
    for g in usable:
        a, b = log_cpm.loc[is_a, g].to_numpy(), log_cpm.loc[~is_a, g].to_numpy()
        p = _welch_p(a, b)
        rows.append(
            {
                "gene": g,
                "arm": CLOCK_PANEL[g],
                "in_index": g in out["index_genes"],
                "mean_log2cpm_a": float(a.mean()),
                "mean_log2cpm_b": float(b.mean()),
                "diff_a_minus_b": float(a.mean() - b.mean()),
                "p": 1.0 if math.isnan(p) else p,
            }
        )
    per_gene = pd.DataFrame(rows).set_index("gene")
    from recluster.evaluators.de import bh_adjust

    per_gene["padj"] = bh_adjust(per_gene["p"].to_numpy())
    n_de = int((per_gene["padj"] < thr.alpha).sum())
    different = [str(g) for g in per_gene.index[per_gene["padj"] < thr.alpha]]

    # phase-structured index across units
    idx = index.to_numpy()

    def stat(v: np.ndarray, mask: np.ndarray) -> float:
        return abs(float(v[mask].mean() - v[~mask].mean()))

    perm = _permutation_p(idx, is_a, stat)
    d = float(idx[is_a].mean() - idx[~is_a].mean())
    pooled = (
        math.sqrt((idx[is_a].var(ddof=1) + idx[~is_a].var(ddof=1)) / 2)
        if n_a > 1 and n_b > 1
        else 0
    )
    can_reject = bool(perm["min_attainable_p"] < thr.alpha)
    out.update(
        {
            "assessable": True,
            "per_gene": per_gene.reset_index().to_dict(orient="records"),
            "n_panel_genes_different": n_de,
            "genes_different": different,
            "index_difference_a_minus_b": d,
            "index_cohens_d": d / pooled if pooled > 0 else None,
            "index_p_permutation": perm["p"],
            "index_n_permutations": perm["n_permutations"],
            "index_exact": perm["exact"],
            "min_attainable_p": perm["min_attainable_p"],
            "test_can_reject": can_reject,
            "index_welch_p": _welch_p(idx[is_a], idx[~is_a]),
            "index_by_sample": {s: float(v) for s, v in index.items()},
            "coherent_shift": bool(can_reject and perm["p"] < thr.alpha),
            "any_gene_differs": bool(n_de > 0),
        }
    )
    return out


# -------------------------------------------------------------- rhythmic-gene list


RHYTHMIC_MIN_TISSUES = 5
RHYTHMIC_MAX_P = 0.05


def rhythmic_overlap(
    genes: list[str], tissue: str | None, table: pd.DataFrame | None
) -> dict[str, Any]:
    """Claimed genes with documented day-night variation in the public table.

    The public table (Wucher et al. 2023, S8) is gene-level: for each gene it gives the number
    of GTEx tissues in which it is higher by day or by night and a binomial P value, but not the
    per-tissue membership, so the check cannot be tissue-specific. A gene is flagged when it is
    day-night in at least ``RHYTHMIC_MIN_TISSUES`` tissues with binomial P below
    ``RHYTHMIC_MAX_P`` (a consistent direction across tissues).
    """
    if table is None or table.empty:
        return {"available": False, "note": "no public day-night gene list available; step skipped"}
    sel = table[
        (table["total_tissues"] >= RHYTHMIC_MIN_TISSUES)
        & (table["p_binomial_all_tissues"] < RHYTHMIC_MAX_P)
    ]
    hits = sel[sel["gene"].isin(genes)]
    flagged = sorted(hits["gene"].astype(str))
    note = (
        f"genes day-night in >= {RHYTHMIC_MIN_TISSUES} of 46 GTEx tissues with a consistent "
        f"direction (binomial P < {RHYTHMIC_MAX_P}); the public table has no per-tissue lists, "
        "so the check is not tissue-specific"
        + (f"; requested tissue {tissue!r} could not be applied" if tissue else "")
    )
    return {
        "available": True,
        "tissue": tissue,
        "flagged": flagged,
        "details": {
            str(r.gene): {
                "day_tissues": int(r.day_tissues),
                "night_tissues": int(r.night_tissues),
                "p": float(r.p_binomial_all_tissues),
            }
            for r in hits.itertuples()
        },
        "note": note,
    }


# ------------------------------------------------------------------------ verdict


@dataclass
class TimeCheck:
    mode: str  # "A" | "B" | "none"
    label: str  # NO TIME CONFOUND DETECTED | POSSIBLE TIME CONFOUND | NOT ASSESSABLE
    reasons: list[str] = field(default_factory=list)
    headline: str = ""  # one sentence for the verdict summary
    balance: dict[str, Any] | None = None
    adjustment: dict[str, Any] | None = None
    panel: dict[str, Any] | None = None
    rhythmic: dict[str, Any] | None = None
    effect_removed: bool = False
    effect_shrunk: bool = False
    disclaimer: str = (
        "This check screens for a collection-time confound at the sample level. It does not "
        "estimate circadian phase from single-cell data."
    )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _direction_ok(value: float | None, direction: str, *, positive: str, negative: str) -> bool:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    if direction == positive:
        return value > 0
    if direction == negative:
        return value < 0
    return True


def _combine(values: list[bool | None], require: str) -> bool | None:
    """all/any over three-valued outcomes; None means 'could not be determined'."""
    if require == "all":
        if any(v is False for v in values):
            return False
        return None if any(v is None for v in values) else bool(values)
    if any(v is True for v in values):
        return True
    return None if any(v is None for v in values) else False


def summarise_adjustment(
    df: pd.DataFrame,
    require: str,
    thr: Thresholds,
    *,
    direction: str = "either",
    info: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Per-gene survival of the claim after adjustment, and the headline flags.

    A gene whose adjusted p or estimate is undefined gets ``survives = None`` and the combined
    outcome is None rather than False: an undefined test is not a vanished effect.
    """
    from recluster.evaluators.de import bh_adjust

    info = info or {}
    out: dict[str, Any] = {
        "design": info.get("design", "~condition"),
        "estimable": bool(info.get("estimable", True)),
        "reason": info.get("reason"),
        "per_gene": [],
        "held_unadjusted": None,
        "survives_adjustment": None,
        "shrunk": False,
        "flipped": False,
    }
    if not out["estimable"]:
        return out
    df = df.copy()
    df["padj_unadjusted"] = bh_adjust(df["p_unadjusted"].to_numpy())
    df["padj_adjusted"] = bh_adjust(df["p_adjusted"].to_numpy())
    rows = []
    for g, r in df.iterrows():
        lu, la = _f(r["lfc_unadjusted"]), _f(r["lfc_adjusted"])
        pu, pa = _f(r["padj_unadjusted"]), _f(r["padj_adjusted"])
        shrink = None if lu is None or la is None or lu == 0 else float(1 - abs(la) / abs(lu))
        held = bool(
            pu is not None
            and pu < thr.alpha
            and _direction_ok(lu, direction, positive="up_in_a", negative="down_in_a")
        )
        survives: bool | None
        if pa is None or la is None or lu is None:
            survives = None
        else:
            survives = bool(
                pa < thr.alpha
                and np.sign(lu) == np.sign(la)
                and _direction_ok(la, direction, positive="up_in_a", negative="down_in_a")
            )
        rows.append(
            {
                "gene": g,
                "lfc_unadjusted": lu,
                "padj_unadjusted": pu,
                "lfc_adjusted": la,
                "padj_adjusted": pa,
                "shrinkage": shrink,
                "direction_flipped": bool(
                    lu is not None and la is not None and np.sign(lu) != np.sign(la) and la != 0
                ),
                "survives": survives,
                "held_unadjusted": held,
            }
        )
    out["per_gene"] = rows
    out["held_unadjusted"] = bool(_combine([r["held_unadjusted"] for r in rows], require))
    out["survives_adjustment"] = _combine([r["survives"] for r in rows], require)
    out["shrunk"] = any(
        r["shrinkage"] is not None and r["shrinkage"] >= thr.time_shrink_flag for r in rows
    )
    out["flipped"] = any(r["direction_flipped"] for r in rows)
    return out


def _f(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def _balance_sentence(balance: dict[str, Any], thr: Thresholds) -> str:
    s = (
        "circular mean collection times differ by "
        f"{balance['circular_mean_difference_hours']:+.1f} h (A minus B); "
        f"permutation p = {balance['p_permutation']:.3g}"
    )
    if not balance.get("test_can_reject", True):
        s += (
            f" (with {balance.get('n_units_a')} vs {balance.get('n_units_b')} units the exact "
            f"test cannot go below {balance.get('min_attainable_p', 0):.3g})"
        )
    return s


def decide(
    mode: str,
    balance: dict[str, Any] | None,
    adjustment: dict[str, Any] | None,
    panel: dict[str, Any] | None,
    thr: Thresholds,
    *,
    held_unadjusted: bool | None = None,
) -> TimeCheck:
    """Apply the documented rules to the computed pieces.

    ``held_unadjusted`` is the baseline run's own outcome when known; it takes precedence over
    the refit's unadjusted result so that the time check never contradicts the grid.
    """
    tc = TimeCheck(
        mode=mode, label=NOT_ASSESSABLE, balance=balance, adjustment=adjustment, panel=panel
    )
    if mode == "A":
        _decide_mode_a(tc, balance, adjustment, thr, held_unadjusted)
    elif mode == "B":
        _decide_mode_b(tc, panel, thr)
    else:
        tc.reasons.append("claim has no condition comparison")
        tc.headline = "not applicable"
    return tc


def _decide_mode_a(
    tc: TimeCheck,
    balance: dict[str, Any] | None,
    adjustment: dict[str, Any] | None,
    thr: Thresholds,
    held_unadjusted: bool | None,
) -> None:
    if not balance or not balance.get("assessable"):
        why = balance.get("reason", "balance not assessable") if balance else "no times"
        tc.reasons.append(why)
        tc.headline = f"not assessable ({why})"
        return
    imbalanced = balance.get("imbalanced")
    if imbalanced is None:  # pieces assembled by hand: fall back to the two rules
        imbalanced = bool(
            balance["p_permutation"] < thr.alpha
            or abs(balance.get("circular_mean_difference_hours", 0.0)) >= thr.time_imbalance_hours
        )
    tc.reasons.append(_balance_sentence(balance, thr))
    if not imbalanced:
        tc.label = NO_CONFOUND
        tc.reasons.append(
            "collection time is balanced between the groups (difference below "
            f"{thr.time_imbalance_hours:g} h and permutation p at or above {thr.alpha:g})"
        )
        tc.headline = "collection time is balanced between the groups"
        return
    basis = balance.get("imbalance_basis") or "permutation test"
    tc.reasons.append(f"collection time is associated with the groups ({basis})")
    if not adjustment:
        tc.label = POSSIBLE
        tc.reasons.append(
            "the claim could not be refitted with time covariates (no population cells)"
        )
        tc.headline = (
            "collection time is associated with the groups; the claim could not be refitted"
        )
        return
    if not adjustment.get("estimable", True):
        tc.label = POSSIBLE
        why = adjustment.get("reason") or "the adjusted model is not estimable"
        tc.reasons.append(f"the effect cannot be separated from collection time: {why}")
        tc.headline = "collection time is associated with the groups and inseparable from condition"
        return
    held = adjustment.get("held_unadjusted") if held_unadjusted is None else held_unadjusted
    if held is False:
        tc.label = POSSIBLE
        tc.reasons.append(
            "the claim did not hold before adjustment, so the adjustment is uninformative; the "
            "design remains confounded with time"
        )
        tc.headline = "collection time is associated with the groups; the claim did not hold anyway"
        return
    survives = adjustment.get("survives_adjustment")
    if survives is None:
        tc.label = POSSIBLE
        tc.reasons.append(
            "the adjusted estimate is undefined for a claimed gene (for example no counts), so "
            "the effect cannot be separated from collection time"
        )
        tc.headline = (
            "collection time is associated with the groups; the adjusted estimate is undefined"
        )
        return
    if survives is False:
        if adjustment.get("shrunk") or adjustment.get("flipped"):
            tc.label, tc.effect_removed = POSSIBLE, True
            change = (
                "its direction flips"
                if adjustment.get("flipped")
                else f"its estimate shrinks by at least {thr.time_shrink_flag:.0%}"
            )
            tc.reasons.append(
                "after adjusting for cos/sin of collection time the effect is no longer "
                f"significant and {change}: the claimed effect is explained by time"
            )
            tc.headline = "the effect is explained by collection time"
            return
        tc.label = POSSIBLE
        tc.reasons.append(
            "after adjusting for cos/sin of collection time the effect is no longer significant "
            "although its estimate changes little (the adjusted model spends degrees of freedom "
            "on time); time and condition cannot be separated with this design"
        )
        tc.headline = "collection time is associated with the groups and inseparable from condition"
        return
    if adjustment.get("shrunk"):
        tc.label, tc.effect_shrunk = POSSIBLE, True
        tc.reasons.append(
            f"the effect survives adjustment but its estimate shrinks by at least "
            f"{thr.time_shrink_flag:.0%}"
        )
        tc.headline = "the effect shrinks substantially after adjusting for collection time"
        return
    tc.label = NO_CONFOUND
    tc.reasons.append(
        "collection time is associated with the groups, but the effect survives adjustment "
        "with little change"
    )
    tc.headline = "the effect survives adjustment for collection time with little change"


def _decide_mode_b(tc: TimeCheck, panel: dict[str, Any] | None, thr: Thresholds) -> None:
    if not panel or not panel.get("assessable"):
        why = panel.get("reason", "panel not assessable") if panel else "no panel"
        tc.reasons.append(why)
        tc.headline = f"not assessable ({why})"
        return
    n_use = len(panel.get("genes_usable") or [])
    n_de = panel.get("n_panel_genes_different", 0)
    genes_de = ", ".join(panel.get("genes_different") or [])
    if not panel.get("test_can_reject", True):
        tc.label = NOT_ASSESSABLE
        welch = panel.get("index_welch_p")
        tc.reasons.append(
            f"{panel.get('n_samples_a')} vs {panel.get('n_samples_b')} units: the exact "
            "permutation test on the clock-state index cannot reach "
            f"{thr.alpha:g} (minimum p = {panel.get('min_attainable_p', 0):.3g}); index "
            f"difference {panel.get('index_difference_a_minus_b', 0):+.2f}"
            + (f", Welch p = {welch:.3g}" if welch is not None and not math.isnan(welch) else "")
            + " reported for information"
        )
        if n_de:
            tc.reasons.append(f"{n_de} of {n_use} clock genes differ gene by gene ({genes_de})")
        tc.headline = "too few units for the sample-level clock-gene test"
        return
    if panel["coherent_shift"]:
        tc.label = POSSIBLE
        tc.reasons.append(
            f"{WARNING_TEXT} (clock-state index differs, permutation p = "
            f"{panel['index_p_permutation']:.3g}; {n_de} of {n_use} panel genes differ "
            "gene by gene" + (f": {genes_de}" if genes_de else "") + ")"
        )
        tc.headline = "clock-gene state differs between the groups in the pattern of a time shift"
        return
    tc.label = NO_CONFOUND
    if panel["any_gene_differs"]:
        tc.reasons.append(
            f"{n_de} of {n_use} clock genes differ between groups gene by gene ({genes_de}), but "
            "not in the antiphase pattern a collection-time difference produces (index p = "
            f"{panel['index_p_permutation']:.3g}); consistent with a direct effect on clock "
            "genes, which this screen cannot distinguish from time"
        )
        tc.headline = (
            f"{n_de} clock gene(s) differ between the groups, not in the pattern of a time shift"
        )
    else:
        tc.reasons.append("no clock-gene difference between groups at the sample level")
        tc.headline = "no clock-gene difference between the groups"


# ------------------------------------------------------------- claim-level entry


def _composition_adjustment_summary(
    adj: dict[str, Any], direction: str, thr: Thresholds
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "design": adj.get("design", "unpaired"),
        "estimable": bool(adj.get("estimable", True)),
        "reason": adj.get("reason"),
        "per_gene": [],
        "effect_unadjusted": _f(adj.get("effect_unadjusted")),
        "p_unadjusted": _f(adj.get("p_unadjusted")),
        "effect_adjusted": _f(adj.get("effect_adjusted")),
        "p_adjusted": _f(adj.get("p_adjusted")),
        "shrinkage": None,
        "held_unadjusted": None,
        "survives_adjustment": None,
        "shrunk": False,
        "flipped": False,
    }
    if not out["estimable"]:
        return out
    b0, p0, b1, p1 = (
        out["effect_unadjusted"],
        out["p_unadjusted"],
        out["effect_adjusted"],
        out["p_adjusted"],
    )
    held = bool(
        p0 is not None
        and p0 < thr.alpha
        and _direction_ok(b0, direction, positive="higher_in_a", negative="lower_in_a")
    )
    survives: bool | None
    if p1 is None or b1 is None or b0 is None:
        survives = None
    else:
        survives = bool(
            p1 < thr.alpha
            and np.sign(b1) == np.sign(b0)
            and _direction_ok(b1, direction, positive="higher_in_a", negative="lower_in_a")
        )
    shrink = None if b0 is None or b1 is None or b0 == 0 else float(1 - abs(b1) / abs(b0))
    out.update(
        {
            "shrinkage": shrink,
            "held_unadjusted": held,
            "survives_adjustment": survives,
            "shrunk": bool(shrink is not None and shrink >= thr.time_shrink_flag),
            "flipped": bool(
                b0 is not None and b1 is not None and np.sign(b0) != np.sign(b1) and b1 != 0
            ),
        }
    )
    return out


def load_rhythmic_table() -> pd.DataFrame | None:
    """The committed day-night gene table built by scripts/fetch_daynight_genes.py, or None."""
    from importlib.resources import files

    try:
        text = files("recluster").joinpath("resources", "daynight_genes.csv").read_text()
    except (FileNotFoundError, OSError, TypeError):
        return None
    df = pd.read_csv(io.StringIO(text), comment="#")
    needed = {"gene", "total_tissues", "day_tissues", "night_tissues", "p_binomial_all_tissues"}
    return df if needed <= set(df.columns) else None


def assess_claim(
    spec: Any,
    counts: AnnData,
    sample: pd.Series,
    condition: pd.Series,
    obs: pd.DataFrame,
    cells: pd.Index | None,
    run_cells: pd.Index,
    thr: Thresholds,
    *,
    n_cpus: int,
    rhythmic_table: pd.DataFrame | None = None,
    pair: pd.Series | None = None,
    held_unadjusted: bool | None = None,
) -> TimeCheck:
    """Run the applicable mode for one claim and return the labelled result.

    ``cells`` are the population cells of the baseline run (or the labelled cells when the
    baseline found nothing), ``pair`` the verdict's pair key, and ``held_unadjusted`` the
    baseline run's own outcome for the claim.
    """
    cmp = spec.comparison
    if cmp is None:
        return decide("none", None, None, None, thr)
    cols = spec.columns
    warns: list[str] = []
    balance = adjustment = panel = None
    hours: pd.Series | None = None
    mode = "A" if cols.collection_time else "B"
    if mode == "A":
        hours, warns = sample_times(obs[cols.collection_time], sample)
    unit, unit_col, notes = exchangeable_units(
        obs,
        counts.obs_names,
        sample,
        condition,
        cmp.group_a,
        cmp.group_b,
        cols.sample_key,
        cols.condition,
        hours_by_sample=hours,
    )
    group_by_unit = condition.loc[unit.index].astype(str).groupby(unit.to_numpy()).first()
    if mode == "A" and hours is not None:
        unit_of_sample = _unit_of_sample(unit, sample)
        hours_by_unit = (
            hours.groupby(unit_of_sample.reindex(hours.index).to_numpy())
            .agg(_circular_mean_of_hours)
            .astype(float)
        )
        balance = balance_check(
            hours_by_unit, group_by_unit, cmp.group_a, cmp.group_b, thr, unit_column=unit_col
        )
        if balance.get("imbalanced") and cells is not None and len(cells):
            if spec.claim_type == "differential_expression":
                df, info = adjusted_pseudobulk_de(
                    counts,
                    cells,
                    sample,
                    condition,
                    hours,
                    cmp,
                    spec.genes,
                    thr,
                    n_cpus=n_cpus,
                    pair=pair,
                )
                adjustment = summarise_adjustment(
                    df, spec.require, thr, direction=cmp.direction, info=info
                )
            elif spec.claim_type == "composition_shift":
                from recluster.evaluators.composition import per_sample_proportions

                table = per_sample_proportions(cells, run_cells, sample, condition, pair)
                adjustment = _composition_adjustment_summary(
                    adjusted_composition(table, hours, cmp.group_a, cmp.group_b),
                    cmp.direction,
                    thr,
                )
    else:
        panel = panel_check(
            counts, unit, group_by_unit, cmp.group_a, cmp.group_b, thr, unit_column=unit_col
        )
    tc = decide(mode, balance, adjustment, panel, thr, held_unadjusted=held_unadjusted)
    tc.reasons = warns + notes + tc.reasons
    if spec.genes:
        tc.rhythmic = rhythmic_overlap(list(spec.genes), spec.tissue, rhythmic_table)
    return tc
