"""Cluster-level doublet check.

Scrublet scores cells one at a time and misses a block of similar doublets that are each
other's neighbours. This check looks at the population as a whole: is its expression centroid
better explained as a point on the segment between two other clusters' centroids than as a
population of its own, and is its library size larger than either parent's? Both are expected
of droplets that contain two cells (one from each parent) and of nothing else that normally
forms a cluster. The idea follows the cluster-based doublet calls in scDblFinder.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData


@dataclass
class DoubletClusterResult:
    best_single: str | None
    single_residual: float | None
    best_pair: list[str] | None
    pair_weight: float | None
    pair_residual: float | None
    residual_ratio: float | None
    libsize_ratio: float | None
    libsize_ratio_vs_max: float | None
    flagged: bool
    note: str = ""

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def cluster_profiles(
    counts: AnnData, clusters: pd.Series, *, n_hvg: int = 2000
) -> tuple[pd.DataFrame, pd.Series]:
    """Per-cluster centroids over the top variable genes (CP10k log1p) and median library size."""
    cells = clusters.index
    a = counts[counts.obs_names.get_indexer(cells)].copy()
    libsize = pd.Series(np.asarray(a.X.sum(axis=1)).ravel(), index=cells)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars - 1), flavor="seurat")
    hv = a.var["highly_variable"].to_numpy()
    x = a.X[:, hv]
    x = sp.csr_matrix(x) if sp.issparse(x) else sp.csr_matrix(np.asarray(x))
    labels = clusters.to_numpy().astype(str)
    uniq, inverse = np.unique(labels, return_inverse=True)
    indicator = sp.csr_matrix(
        (np.ones(len(labels)), (inverse, np.arange(len(labels)))), shape=(len(uniq), len(labels))
    )
    sums = indicator @ x
    sums = sums.toarray() if sp.issparse(sums) else np.asarray(sums)
    sizes = np.bincount(inverse, minlength=len(uniq)).astype(float)
    centroids = pd.DataFrame(sums / sizes[:, None], index=uniq)
    med = pd.Series(libsize.to_numpy(), index=cells).groupby(labels).median()
    return centroids, med.reindex(uniq)


def _segment_fit(t: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Least-squares projection of t onto the segment a-b; returns (weight on a, residual)."""
    d = a - b
    denom = float(d @ d)
    if denom == 0:
        return 1.0, float(np.linalg.norm(t - a))
    w = float(np.clip(((t - b) @ d) / denom, 0.0, 1.0))
    return w, float(np.linalg.norm(t - (w * a + (1 - w) * b)))


def mixture_fit(
    target: np.ndarray,
    target_libsize: float,
    others: pd.DataFrame,
    other_libsize: pd.Series,
) -> DoubletClusterResult:
    """Compare the best single-cluster explanation of ``target`` with the best two-cluster mix."""
    names = list(others.index)
    if len(names) < 2:
        return DoubletClusterResult(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            False,
            "fewer than 2 other clusters",
        )
    mat = others.to_numpy()
    singles = np.linalg.norm(mat - target[None, :], axis=1)
    i_single = int(np.argmin(singles))
    best_pair: tuple[int, int] | None = None
    best_w, best_res = 0.0, float("inf")
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            w, res = _segment_fit(target, mat[i], mat[j])
            if res < best_res:
                best_pair, best_w, best_res = (i, j), w, res
    assert best_pair is not None
    i, j = best_pair
    single_res = float(singles[i_single])
    ratio = best_res / single_res if single_res > 0 else float("nan")
    lib_mix = best_w * float(other_libsize.iloc[i]) + (1 - best_w) * float(other_libsize.iloc[j])
    lib_max = max(float(other_libsize.iloc[i]), float(other_libsize.iloc[j]))
    return DoubletClusterResult(
        best_single=names[i_single],
        single_residual=single_res,
        best_pair=[names[i], names[j]],
        pair_weight=best_w,
        pair_residual=best_res,
        residual_ratio=ratio,
        libsize_ratio=target_libsize / lib_mix if lib_mix > 0 else None,
        libsize_ratio_vs_max=target_libsize / lib_max if lib_max > 0 else None,
        flagged=False,
    )


def doublet_cluster_check(
    counts: AnnData,
    clusters: pd.Series,
    population_clusters: list[str],
    *,
    residual_ratio_max: float,
    libsize_ratio_min: float,
) -> DoubletClusterResult:
    """Run the mixture test for the matched population (a union of clusters) against the rest."""
    if not population_clusters:
        return DoubletClusterResult(
            None, None, None, None, None, None, None, None, False, "no population"
        )
    centroids, libsize = cluster_profiles(counts, clusters)
    labels = clusters.to_numpy().astype(str)
    in_pop = np.isin(labels, population_clusters)
    pop_cells = clusters.index[in_pop]
    # centroid of the union, weighted by cluster size
    sizes = pd.Series(labels[in_pop]).value_counts()
    pop_centroid = np.average(
        centroids.loc[sizes.index].to_numpy(), axis=0, weights=sizes.to_numpy()
    )
    pop_lib = float(
        np.median(np.asarray(counts[counts.obs_names.get_indexer(pop_cells)].X.sum(axis=1)).ravel())
    )
    others = centroids.drop(index=population_clusters, errors="ignore")
    result = mixture_fit(pop_centroid, pop_lib, others, libsize.reindex(others.index))
    if result.residual_ratio is not None and result.libsize_ratio is not None:
        result.flagged = bool(
            result.residual_ratio <= residual_ratio_max
            and result.libsize_ratio >= libsize_ratio_min
        )
    return result
