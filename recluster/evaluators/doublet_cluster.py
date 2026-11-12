"""Cluster-level doublet check.

Scrublet scores cells one at a time and misses a block of similar doublets that are each
other's neighbours. This check looks at the population as a whole: is its expression centroid
better explained as a point on the segment between two other clusters' centroids than as a
population of its own, and is its library size larger than either parent's? Both are expected
of droplets that contain two cells (one from each parent) and of nothing else that normally
forms a cluster. The idea follows the cluster-based doublet calls in scDblFinder. It is
reported for the matched population and does not change the verdict.

``doublet_population_test`` is the removal step the pipeline uses: every cluster against
simulated doublets of other cell types, with a cross-validated classifier (PLAN.md, v0.6).
"""

from __future__ import annotations

import warnings
import zlib
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


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


# --------------------------------------------------------------------------------------------
# Population-level two-sample test (PLAN.md "v0.6 round 2"), used by the pipeline to remove
# clusters made of heterotypic doublets when doublet removal is on.

_N_HVG = 2000
_N_PCS = 20
_N_SHORTLIST = 3
_N_REP = 5
_N_FOLDS = 5
_MAX_C = 500
_GROUP_CUT = 1.0


@dataclass
class DoubletPopulationRecord:
    """Test result for one cluster.

    ``statistic`` is None when the other clusters form fewer than two groups (no parent pair).
    ``best_pair`` gives the member clusters of the two parent groups behind the statistic.
    """

    cluster: str
    size: int
    statistic: float | None
    best_pair: tuple[list[str], list[str]] | None
    flagged: bool


@dataclass
class DoubletPopulationResult:
    records: list[DoubletPopulationRecord]
    flagged: list[str]
    auroc_max: float

    def statistics(self) -> dict[str, float | None]:
        return {r.cluster: r.statistic for r in self.records}


def _key(*parts: object) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode())


def _rng(seed: int, *parts: object) -> np.random.Generator:
    return np.random.default_rng([seed, _key(*parts)])


def _cluster_order(name: str) -> tuple[int, int, str]:
    """Numeric cluster ids in numeric order, then the rest alphabetically."""
    return (0, int(name), "") if name.isdigit() else (1, 0, name)


def _lognorm(raw: sp.csr_matrix, hvg: np.ndarray) -> np.ndarray:
    """log1p of counts per 10k over the given genes (library size over all genes)."""
    tot = np.asarray(raw.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    return np.log1p(raw[:, hvg].toarray() / tot[:, None] * 1e4)


def _size_features(raw: sp.csr_matrix) -> np.ndarray:
    """log10 library size and log10 genes detected per cell."""
    lib = np.log10(np.maximum(np.asarray(raw.sum(axis=1)).ravel(), 1))
    n_genes = np.log10(np.maximum(np.diff(raw.indptr), 1))
    return np.column_stack([lib, n_genes])


def _oof_auroc(
    xa: np.ndarray, sa: np.ndarray, xb: np.ndarray, sb: np.ndarray, pca_seed: int, cv_seed: int
) -> float:
    """Out-of-fold AUROC of a logistic regression telling set a from set b.

    Features: PCs of the log expression fitted on both sets, plus the two size features.
    """
    x = np.vstack([xa, xb])
    y = np.r_[np.zeros(len(xa)), np.ones(len(xb))]
    n_pcs = min(_N_PCS, x.shape[0] - 1, x.shape[1])
    pcs = PCA(n_components=n_pcs, svd_solver="randomized", random_state=pca_seed).fit_transform(x)
    feats = np.column_stack([pcs, np.vstack([sa, sb])])
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000))
    cv = StratifiedKFold(_N_FOLDS, shuffle=True, random_state=cv_seed)
    scores = cross_val_predict(clf, feats, y, cv=cv, method="decision_function")
    return float(roc_auc_score(y, scores))


def _coarse_groups(zcent: np.ndarray) -> np.ndarray:
    """Average-linkage groups of z-scored centroids, 1 - Pearson distance cut at 1.0."""
    if len(zcent) == 1:
        return np.array([1])
    lk = linkage(pdist(zcent, "correlation"), "average")
    return fcluster(lk, t=_GROUP_CUT, criterion="distance")


def doublet_population_test(
    counts: AnnData,
    clusters: pd.Series,
    *,
    seed: int,
    auroc_max: float,
    min_cluster_size: int = 10,
) -> DoubletPopulationResult:
    """Flag clusters that cannot be told apart from simulated doublets of two other cell types.

    ``counts`` holds raw counts (all genes); ``clusters`` maps barcodes of one run's kept cells to
    cluster ids. For every cluster C with at least ``min_cluster_size`` cells:

    * the other such clusters are joined into cell-type groups by average-linkage clustering of
      their centroids (per-gene z-scored log expression clipped at +-10, 1 - Pearson, cut 1.0);
    * the 3 pairs of groups whose centroid segment lies closest to C's centroid are the
      candidate parents;
    * for each pair, ``min(max(|C|, 100), 500)`` doublets are simulated by summing the raw counts
      of a random cell of each group, and a cross-validated classifier separates C (subsampled to
      500 cells) from them; the out-of-fold AUROC is averaged over 5 repetitions;
    * the statistic is the minimum over the pairs; C is flagged when it is below ``auroc_max``.

    Expression is log1p of counts per 10k over 2000 Seurat-flavour variable genes chosen on these
    cells. Everything random is seeded from ``seed`` and the cluster and pair names.
    """
    labels = clusters.astype(str)
    rows = counts.obs_names.get_indexer(labels.index)
    if (rows < 0).any():
        raise KeyError(f"{int((rows < 0).sum())} clustered cells are missing from the counts")
    raw = sp.csr_matrix(counts.X[rows], dtype=np.float32)
    norm = AnnData(raw.copy())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(norm, target_sum=1e4)
        sc.pp.log1p(norm)
        sc.pp.highly_variable_genes(norm, n_top_genes=min(_N_HVG, norm.n_vars), flavor="seurat")
    hvg = np.flatnonzero(norm.var["highly_variable"].to_numpy())
    x_all = sp.csr_matrix(norm.X[:, hvg], dtype=np.float64)
    del norm
    mu = np.asarray(x_all.mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(np.asarray(x_all.multiply(x_all).mean(axis=0)).ravel() - mu**2, 0))
    sd[sd == 0] = 1
    s_raw = _size_features(raw)
    s_mu, s_sd = s_raw.mean(axis=0), s_raw.std(axis=0)
    s_sd[s_sd == 0] = 1
    s_all = (s_raw - s_mu) / s_sd

    lab = labels.to_numpy()
    sizes = pd.Series(lab).value_counts()
    big = sorted((str(k) for k in sizes[sizes >= min_cluster_size].index), key=_cluster_order)
    idx = {k: np.flatnonzero(lab == k) for k in big}
    cent: dict[str, np.ndarray] = {}
    zcent: dict[str, np.ndarray] = {}
    for k in big:
        xk = x_all[idx[k]].toarray()
        cent[k] = xk.mean(axis=0)
        zcent[k] = np.clip((xk - mu) / sd, -10, 10).mean(axis=0)

    def simulate(
        rng: np.random.Generator, a: np.ndarray, b: np.ndarray, n: int
    ) -> tuple[np.ndarray, np.ndarray]:
        dbl = sp.csr_matrix(
            raw[rng.choice(a, n, replace=True)] + raw[rng.choice(b, n, replace=True)]
        )
        return _lognorm(dbl, hvg), (_size_features(dbl) - s_mu) / s_sd

    records: list[DoubletPopulationRecord] = []
    for c in big:
        n_c = len(idx[c])
        others = [k for k in big if k != c]
        pairs: list[tuple[float, int, int]] = []
        members: dict[int, list[str]] = {}
        if len(others) >= 2:
            g = _coarse_groups(np.array([zcent[k] for k in others]))
            gids = sorted({int(v) for v in g})
            members = {gi: [k for k, gg in zip(others, g, strict=True) if gg == gi] for gi in gids}
            gcent = {
                gi: np.average(
                    [cent[k] for k in members[gi]],
                    axis=0,
                    weights=[len(idx[k]) for k in members[gi]],
                )
                for gi in gids
            }
            for i, ga in enumerate(gids):
                for gb in gids[i + 1 :]:
                    _, res = _segment_fit(cent[c], gcent[ga], gcent[gb])
                    pairs.append((res, ga, gb))
        if not pairs:
            records.append(DoubletPopulationRecord(c, n_c, None, None, False))
            continue
        pairs.sort(key=lambda p: p[0])
        n_sim = min(max(n_c, 100), 500)
        if n_c <= _MAX_C:
            csub = idx[c]
        else:
            csub = np.sort(_rng(seed, c, "csub").choice(idx[c], _MAX_C, replace=False))
        xc, sc_ = x_all[csub].toarray(), s_all[csub]
        best: tuple[float, int, int] | None = None
        for _, ga, gb in pairs[:_N_SHORTLIST]:
            gidx_a = np.concatenate([idx[k] for k in members[ga]])
            gidx_b = np.concatenate([idx[k] for k in members[gb]])
            name = ",".join(members[ga]) + "+" + ",".join(members[gb])
            aucs = []
            for rep in range(_N_REP):
                xs, ss = simulate(_rng(seed, c, name, rep, "sim"), gidx_a, gidx_b, n_sim)
                cv_seed = _key(seed, c, name, rep, "cv") % (2**31)
                aucs.append(_oof_auroc(xc, sc_, xs, ss, pca_seed=seed, cv_seed=cv_seed))
            stat = float(np.mean(aucs))
            if best is None or stat < best[0]:
                best = (stat, ga, gb)
        assert best is not None
        stat, ga, gb = best
        records.append(
            DoubletPopulationRecord(
                c, n_c, stat, (list(members[ga]), list(members[gb])), stat < auroc_max
            )
        )
    return DoubletPopulationResult(
        records=records, flagged=[r.cluster for r in records if r.flagged], auroc_max=auroc_max
    )
