"""differential_expression: cell-level Wilcoxon vs replicate-level pseudobulk (pydeseq2).

Cell-level tests treat every cell as an independent observation and are the classic source of
false positives in multi-sample single-cell DE (Squair et al. 2021, Nat Commun 12:5692). The
verdict uses the pseudobulk test, blocked on the pair (donor) when the design is paired. The
cell-level result is reported for comparison and to flag pseudoreplication.
"""

from __future__ import annotations

import contextlib
import io
import math
import warnings

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from recluster.config import Thresholds
from recluster.evaluators.pseudobulk import pseudobulk_counts
from recluster.spec import Comparison


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg over the non-NaN entries."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order] * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def cell_level_de(
    lognorm: pd.DataFrame, cells_a: pd.Index, cells_b: pd.Index, genes: list[str]
) -> pd.DataFrame:
    """Two-sided Wilcoxon rank-sum per gene on CP10k log1p values; log2FC of means."""
    xa, xb = lognorm.loc[cells_a, genes], lognorm.loc[cells_b, genes]
    rows = []
    for g in genes:
        a, b = xa[g].to_numpy(), xb[g].to_numpy()
        constant = len(a) and len(b) and np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]
        if len(a) < 2 or len(b) < 2 or constant:
            p = math.nan
        else:
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        ma, mb = (
            float(np.expm1(a).mean()) if len(a) else 0.0,
            float(np.expm1(b).mean()) if len(b) else 0.0,
        )
        rows.append(
            {
                "gene": g,
                "cell_log2fc": math.log2((ma + 1.0) / (mb + 1.0)),
                "cell_p": p,
                "pct_expressing_a": float((a > 0).mean()) if len(a) else math.nan,
                "pct_expressing_b": float((b > 0).mean()) if len(b) else math.nan,
            }
        )
    df = pd.DataFrame(rows).set_index("gene")
    df["cell_padj"] = bh_adjust(df["cell_p"].to_numpy())
    return df


def pseudobulk_de(
    counts: AnnData,
    cells: pd.Index,
    sample: pd.Series,
    condition: pd.Series,
    comparison: Comparison,
    genes: list[str],
    thr: Thresholds,
    *,
    n_cpus: int,
    pair: pd.Series | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """DESeq2 (pydeseq2) on per-sample summed counts within the population.

    With a pair key (every pair level in both groups) the design is ``~ pair + condition``,
    which tests the within-pair contrast; otherwise ``~ condition``. ``info`` records the
    design used (``design``, ``design_formula``) and, when a pair key was given but not every
    pair had a pseudobulk in both groups, why it fell back (``design_note``). The formula is
    None when DESeq2 was not run.
    """
    pb, n_cells = pseudobulk_counts(counts, cells, sample, min_cells=thr.min_cells_per_sample)
    by_sample = sample.loc[cells].astype(str).to_numpy()
    cond = condition.loc[cells].astype(str).groupby(by_sample).first().reindex(pb.index)
    keep = cond.isin([comparison.group_a, comparison.group_b])
    pb, cond, n_cells = pb.loc[keep], cond.loc[keep], n_cells.loc[keep]
    n_a, n_b = int((cond == comparison.group_a).sum()), int((cond == comparison.group_b).sum())
    warns: list[str] = []
    info: dict[str, object] = {
        "n_samples_a": n_a,
        "n_samples_b": n_b,
        "cells_per_sample": n_cells.to_dict(),
        "design": "unpaired",
        "design_formula": None,
        "design_note": None,
        "warnings": warns,
    }
    empty = pd.DataFrame(
        {
            "pb_base_mean": math.nan,
            "pb_log2fc": math.nan,
            "pb_p": math.nan,
            "pb_padj_genomewide": math.nan,
        },
        index=pd.Index(genes, name="gene"),
    )
    if n_a < thr.min_samples_present or n_b < thr.min_samples_present:
        warns.append(
            f"population present in {n_a} vs {n_b} pseudobulk replicates "
            f"(need {thr.min_samples_present} per group); DESeq2 not run"
        )
        return empty, info

    totals = pb.sum(axis=0)
    gene_keep = (totals >= thr.pseudobulk_min_total_counts) | pb.columns.isin(genes)
    pb = pb.loc[:, gene_keep]
    info["n_genes_tested"] = int(pb.shape[1])

    meta = pd.DataFrame(
        {"condition": np.where(cond == comparison.group_a, "grpA", "grpB")}, index=pb.index
    )
    design = "~condition"
    if pair is not None:
        # the pair key covers the two compared groups only; cells of other conditions have no
        # pair id, and their samples were dropped above
        by_cell = pair.astype(str).reindex(cells)
        pair_of_sample = by_cell.groupby(by_sample).first().reindex(pb.index)
        levels = {p: f"p{i}" for i, p in enumerate(sorted(pair_of_sample.dropna().unique()))}
        meta["pair"] = pair_of_sample.map(levels).to_numpy()
        complete = meta.groupby("pair")["condition"].nunique()
        if meta["pair"].notna().all() and (complete == 2).all() and len(complete) >= 2:
            design = "~pair + condition"
            info["design"] = "paired"
        else:
            note = "not every pair has a pseudobulk in both groups; unpaired design used"
            info["design_note"] = note
            warns.append(note)
    info["design_formula"] = design
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=pb, metadata=meta, design=design, quiet=True, n_cpus=n_cpus)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["condition", "grpA", "grpB"], quiet=True, n_cpus=n_cpus)
        st.summary()
        res = st.results_df
    out = pd.DataFrame(
        {
            "pb_base_mean": res["baseMean"],
            "pb_log2fc": res["log2FoldChange"],
            "pb_p": res["pvalue"],
            "pb_padj_genomewide": res["padj"],
        }
    ).reindex(genes)
    out.index.name = "gene"
    return out, info


def evaluate_de(
    *,
    cells: pd.Index,
    counts: AnnData,
    lognorm: pd.DataFrame,
    sample: pd.Series,
    condition: pd.Series,
    comparison: Comparison,
    genes: list[str],
    require: str,
    thr: Thresholds,
    n_cpus: int,
    pair: pd.Series | None = None,
) -> dict[str, object]:
    cond = condition.loc[cells].astype(str)
    cells_a = cells[(cond == comparison.group_a).to_numpy()]
    cells_b = cells[(cond == comparison.group_b).to_numpy()]
    cell_df = cell_level_de(lognorm, cells_a, cells_b, genes)
    pb_df, info = pseudobulk_de(
        counts, cells, sample, condition, comparison, genes, thr, n_cpus=n_cpus, pair=pair
    )
    df = cell_df.join(pb_df)
    df["pb_padj"] = bh_adjust(df["pb_p"].to_numpy())  # multiplicity over the claimed genes
    n_a = int(info["n_samples_a"])  # type: ignore[call-overload]
    n_b = int(info["n_samples_b"])  # type: ignore[call-overload]
    evaluable = n_a >= thr.min_samples_present and n_b >= thr.min_samples_present
    reasons = (
        []
        if evaluable
        else [
            f"population present in only {n_a} vs {n_b} replicate(s) "
            f"(need {thr.min_samples_present} per group); the run cannot test this claim"
        ]
    )

    def direction_ok(lfc: float) -> bool:
        if lfc is None or math.isnan(lfc):
            return False
        if comparison.direction == "up_in_a":
            return lfc > 0
        if comparison.direction == "down_in_a":
            return lfc < 0
        return True

    df["holds_cell"] = [
        bool(not math.isnan(p) and p < thr.alpha and direction_ok(lfc))
        for p, lfc in zip(df["cell_padj"], df["cell_log2fc"], strict=True)
    ]
    df["holds_replicate"] = [
        bool(not math.isnan(p) and p < thr.alpha and direction_ok(lfc))
        for p, lfc in zip(df["pb_padj"], df["pb_log2fc"], strict=True)
    ]
    df["agreement"] = [
        "both" if c and r else "cell_only" if c else "pseudobulk_only" if r else "neither"
        for c, r in zip(df["holds_cell"], df["holds_replicate"], strict=True)
    ]
    combine = all if require == "all" else any
    holds = bool(evaluable and combine(df["holds_replicate"].tolist()))
    holds_cell = bool(combine(df["holds_cell"].tolist()))

    per_gene = df.reset_index().replace({np.nan: None}).to_dict(orient="records")
    first = df.iloc[0]
    return {
        "holds": holds,
        "holds_cell_level": holds_cell,
        "replicate_ok": holds,
        "evaluable": evaluable,
        "reasons": reasons,
        "design": info["design"],
        "design_formula": info["design_formula"],
        "design_note": info["design_note"],
        "n_cells_a": int(len(cells_a)),
        "n_cells_b": int(len(cells_b)),
        "n_samples_a": n_a,
        "n_samples_b": n_b,
        "n_genes_tested_pseudobulk": info.get("n_genes_tested"),
        "cells_per_sample": info["cells_per_sample"],
        "per_gene": per_gene,
        "warnings": list(info.get("warnings") or []),  # type: ignore[call-overload]
        "key_stats": {
            "pb_log2fc": _none_if_nan(first["pb_log2fc"]),
            "pb_padj": _none_if_nan(first["pb_padj"]),
            "cell_padj": _none_if_nan(first["cell_padj"]),
        },
    }


def _none_if_nan(x: object) -> float | None:
    try:
        f = float(x)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f
