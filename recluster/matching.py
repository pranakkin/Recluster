"""Find the claimed population inside one grid run.

Two complementary signals: marker enrichment (AUROC of a per-cell marker score for a cluster
versus everything else) and overlap with a reference cell set (Jaccard over barcodes). Cells
removed by QC in a run count as absent, so QC choices that delete the population show up as a
lower Jaccard rather than being silently ignored.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import rankdata

from recluster.config import Thresholds


def lognorm_matrix(counts: AnnData, genes: list[str]) -> pd.DataFrame:
    """CP10k log1p expression for a few genes, all cells. Independent of any grid choice."""
    if not genes:
        return pd.DataFrame(index=counts.obs_names)
    x = counts.X
    total = np.asarray(x.sum(axis=1)).ravel()
    idx = counts.var_names.get_indexer(genes)
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise KeyError(f"genes not in counts matrix: {missing}")
    sub = x[:, idx]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.where(total[:, None] > 0, sub / total[:, None] * 1e4, 0.0)
    return pd.DataFrame(np.log1p(norm), index=counts.obs_names, columns=genes)


def marker_score(lognorm: pd.DataFrame, markers: list[str], negative: list[str]) -> pd.Series:
    """Mean marker expression minus mean negative-marker expression, z-scored over cells."""
    s = lognorm[markers].mean(axis=1)
    if negative:
        s = s - lognorm[negative].mean(axis=1)
    sd = float(s.std(ddof=0))
    z = (s - s.mean()) / sd if sd > 0 else s * 0.0
    return z.rename("marker_score")


def auroc_by_group(score: np.ndarray, groups: np.ndarray) -> dict[str, float]:
    """AUROC of ``score`` for each group vs the rest, via the Mann-Whitney U identity."""
    ranks = rankdata(score)
    n = len(score)
    out: dict[str, float] = {}
    for g in np.unique(groups):
        mask = groups == g
        n1 = int(mask.sum())
        n0 = n - n1
        if n1 == 0 or n0 == 0:
            out[str(g)] = float("nan")
            continue
        u = ranks[mask].sum() - n1 * (n1 + 1) / 2.0
        out[str(g)] = float(u / (n1 * n0))
    return out


def auroc_two_sets(score: np.ndarray, positive: np.ndarray) -> float:
    """AUROC of ``score`` for the boolean ``positive`` mask vs the rest."""
    n1 = int(positive.sum())
    n0 = len(score) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(score)
    u = ranks[positive].sum() - n1 * (n1 + 1) / 2.0
    return float(u / (n1 * n0))


@dataclass
class Match:
    """Where the claimed population landed in one run."""

    mode: str  # "markers" or "labels"
    clusters: list[str]
    n_cells: int
    found: bool
    best_cluster_auroc: float | None
    marker_auroc: float | None
    jaccard: float | None
    precision: float | None
    recall: float | None
    n_reference: int | None
    n_samples_present: int | None
    max_sample_share: float | None
    cells: pd.Index

    def to_dict(self) -> dict[str, object]:
        d = asdict(self)
        d.pop("cells")
        return d


def match_population(
    clusters: pd.Series,
    *,
    score: pd.Series | None,
    reference: pd.Index | None,
    sample: pd.Series | None,
    min_cells: int,
    thr: Thresholds,
) -> Match:
    """Locate the population in a run by marker enrichment (preferred) or label overlap."""
    run_cells = clusters.index
    labels = clusters.to_numpy().astype(str)
    sizes = pd.Series(labels).value_counts()
    best_auroc: float | None = None

    if score is not None:
        s = score.loc[run_cells].to_numpy()
        aurocs = auroc_by_group(s, labels)
        best_auroc = float(np.nanmax(list(aurocs.values()))) if aurocs else None
        candidates = [
            c
            for c, a in aurocs.items()
            if a >= thr.marker_auroc_min and sizes[c] >= thr.min_cluster_size
        ]
        mode = "markers"
    else:
        if reference is None:
            raise ValueError("label-mode matching needs a reference set")
        in_ref = pd.Series(run_cells.isin(reference), index=run_cells)
        precision_by_cluster = in_ref.groupby(labels).mean()
        candidates = [
            str(c)
            for c, p in precision_by_cluster.items()
            if p >= thr.label_precision_min and sizes[str(c)] >= thr.min_cluster_size
        ]
        mode = "labels"

    in_m = np.isin(labels, candidates)
    cells = run_cells[in_m]
    n_cells = int(len(cells))

    marker_auroc = None
    if score is not None and 0 < n_cells < len(run_cells):
        marker_auroc = auroc_two_sets(score.loc[run_cells].to_numpy(), in_m)

    jaccard = precision = recall = None
    n_reference = None
    if reference is not None:
        n_reference = int(len(reference))
        inter = int(cells.isin(reference).sum())
        union = n_cells + n_reference - inter
        jaccard = inter / union if union else 0.0
        precision = inter / n_cells if n_cells else 0.0
        recall = inter / n_reference if n_reference else 0.0

    n_samples_present = max_share = None
    if sample is not None and n_cells > 0:
        per_sample = sample.loc[cells].value_counts()
        n_samples_present = int((per_sample >= thr.min_cells_per_sample).sum())
        max_share = float(per_sample.max() / n_cells)
    elif sample is not None:
        n_samples_present, max_share = 0, None

    return Match(
        mode=mode,
        clusters=sorted(candidates, key=_cluster_sort_key),
        n_cells=n_cells,
        found=n_cells >= min_cells,
        best_cluster_auroc=best_auroc,
        marker_auroc=marker_auroc,
        jaccard=jaccard,
        precision=precision,
        recall=recall,
        n_reference=n_reference,
        n_samples_present=n_samples_present,
        max_sample_share=max_share,
        cells=cells,
    )


def _cluster_sort_key(c: str) -> tuple[int, str]:
    return (int(c), "") if c.isdigit() else (10**9, c)
