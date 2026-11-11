"""Cluster-level doublet check: segment geometry, the two-cluster edge case, and end to end."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from recluster.config import Thresholds
from recluster.evaluators.doublet_cluster import (
    DoubletClusterResult,
    doublet_cluster_check,
    mixture_fit,
)

THR = Thresholds()


def test_default_thresholds():
    assert THR.doublet_residual_ratio_max == 0.60
    assert THR.doublet_libsize_ratio_min == 1.50


def test_mixture_fit_recovers_a_point_on_the_segment():
    rng = np.random.default_rng(0)
    centroids = pd.DataFrame(rng.uniform(0, 5, size=(3, 20)), index=["A", "B", "C"])
    libsize = pd.Series(1000.0, index=centroids.index)
    target = 0.5 * centroids.loc["A"].to_numpy() + 0.5 * centroids.loc["B"].to_numpy()
    res = mixture_fit(target, 2000.0, centroids, libsize)
    assert sorted(res.best_pair) == ["A", "B"]
    assert res.pair_weight == pytest.approx(0.5, abs=0.05)
    assert res.pair_residual == pytest.approx(0.0, abs=1e-9)
    assert res.residual_ratio < 1e-6
    assert res.best_single in ("A", "B") and res.single_residual > 1.0
    assert res.libsize_ratio == pytest.approx(2.0)
    assert res.libsize_ratio_vs_max == pytest.approx(2.0)
    # the thresholds are applied by doublet_cluster_check, never by the fit itself
    assert res.flagged is False and res.note == ""


def test_mixture_fit_needs_two_other_clusters():
    centroids = pd.DataFrame(np.ones((1, 20)), index=["A"])
    res = mixture_fit(np.zeros(20), 500.0, centroids, pd.Series([1000.0], index=["A"]))
    assert res == DoubletClusterResult(
        None, None, None, None, None, None, None, None, False, "fewer than 2 other clusters"
    )
    assert res.to_dict()["note"] == "fewer than 2 other clusters"


def _tiny(seed: int = 0, *, genuine_fourth: bool = False) -> tuple[AnnData, pd.Series]:
    """Three Poisson populations (60 cells, 200 genes) plus a fourth cluster of 40 cells.

    The populations share a base profile with lognormal per-gene fold changes (sigma 0.5).
    The fourth cluster is either the cell-wise sum of a random A cell and a random B cell (a
    block of doublets) or a population with a profile of its own.
    """
    rng = np.random.default_rng(seed)
    n_genes, n_real, n_extra = 200, 60, 40
    base = rng.gamma(2.0, 5.0, size=n_genes)
    blocks: list[np.ndarray] = []
    labels: list[str] = []
    for name in ("A", "B", "C"):
        mu = base * rng.lognormal(0.0, 0.5, size=n_genes)
        blocks.append(rng.poisson(mu[None, :], size=(n_real, n_genes)))
        labels += [name] * n_real
    if genuine_fourth:
        mu = base * rng.lognormal(0.0, 0.5, size=n_genes)
        blocks.append(rng.poisson(mu[None, :], size=(n_extra, n_genes)))
    else:
        ia = rng.integers(n_real, size=n_extra)
        ib = rng.integers(n_real, size=n_extra)
        blocks.append(blocks[0][ia] + blocks[1][ib])
    labels += ["D"] * n_extra
    x = np.vstack(blocks).astype(np.float32)
    obs = pd.DataFrame({"cluster": labels}, index=[f"c{i}" for i in range(len(labels))])
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    adata = AnnData(X=sp.csr_matrix(x), obs=obs, var=var)
    return adata, adata.obs["cluster"].astype(str)


def _check(adata: AnnData, clusters: pd.Series) -> DoubletClusterResult:
    return doublet_cluster_check(
        adata,
        clusters,
        ["D"],
        residual_ratio_max=THR.doublet_residual_ratio_max,
        libsize_ratio_min=THR.doublet_libsize_ratio_min,
    )


def test_doublet_block_is_flagged_with_its_parents():
    adata, clusters = _tiny()
    res = _check(adata, clusters)
    assert set(res.best_pair) == {"A", "B"}
    assert res.pair_weight == pytest.approx(0.5, abs=0.1)
    assert res.residual_ratio < THR.doublet_residual_ratio_max
    assert res.libsize_ratio > 1.5 and res.libsize_ratio_vs_max > 1.5
    assert res.flagged is True
    assert _check(adata, clusters) == res


def test_genuine_fourth_population_is_not_flagged():
    adata, clusters = _tiny(genuine_fourth=True)
    res = _check(adata, clusters)
    assert res.libsize_ratio == pytest.approx(1.0, abs=0.25)
    assert res.residual_ratio > THR.doublet_residual_ratio_max
    assert res.flagged is False


def test_two_cluster_data_cannot_be_tested():
    adata, clusters = _tiny()
    keep = clusters.isin(["A", "D"]).to_numpy()
    res = _check(adata[keep].copy(), clusters[keep])
    assert res.flagged is False and res.note == "fewer than 2 other clusters"


def test_empty_population_is_not_tested():
    adata, clusters = _tiny()
    res = doublet_cluster_check(
        adata, clusters, [], residual_ratio_max=0.60, libsize_ratio_min=1.50
    )
    assert res.flagged is False and res.note == "no population"
