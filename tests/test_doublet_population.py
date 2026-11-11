"""Population-level doublet test and its removal step in the pipeline.

Synthetic data: three Poisson populations with distinct profiles and per-cell size factors, plus
a block of doublets, each the cell-wise sum of a random A cell and a random B cell. The block
forms its own cluster.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from recluster.config import GridConfig, PipelineParams
from recluster.evaluators.doublet_cluster import (
    DoubletPopulationResult,
    doublet_population_test,
)
from recluster.pipeline import Pipeline

N_REAL, N_DOUBLET = 150, 80
PARAMS = PipelineParams(
    min_genes=10, resolution=0.5, n_pcs=10, batch_correction=False, doublet_removal=True
)


@pytest.fixture(scope="module")
def data() -> AnnData:
    rng = np.random.default_rng(0)
    n_genes = 500
    base = rng.gamma(2.0, 0.5, size=n_genes)
    blocks: list[np.ndarray] = []
    labels: list[str] = []
    for name in ("A", "B", "C"):
        mu = base * rng.lognormal(0.0, 1.0, size=n_genes)
        size = rng.lognormal(0.0, 0.3, size=N_REAL)
        blocks.append(rng.poisson(size[:, None] * mu[None, :]))
        labels += [name] * N_REAL
    ia = rng.integers(N_REAL, size=N_DOUBLET)
    ib = rng.integers(N_REAL, size=N_DOUBLET)
    blocks.append(blocks[0][ia] + blocks[1][ib])
    labels += ["D"] * N_DOUBLET
    x = np.vstack(blocks).astype(np.float32)
    obs = pd.DataFrame({"truth": labels}, index=[f"c{i}" for i in range(len(labels))])
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return AnnData(X=sp.csr_matrix(x), obs=obs, var=var)


@pytest.fixture(scope="module")
def result(data: AnnData) -> DoubletPopulationResult:
    return doublet_population_test(data, data.obs["truth"].astype(str), seed=0, auroc_max=0.80)


def test_doublet_block_is_flagged_and_real_populations_are_not(
    result: DoubletPopulationResult,
) -> None:
    rec = {r.cluster: r for r in result.records}
    assert result.flagged == ["D"]
    assert rec["D"].statistic is not None and rec["D"].statistic < 0.80
    assert rec["D"].size == N_DOUBLET and rec["D"].flagged
    assert rec["D"].best_pair is not None
    assert sorted(rec["D"].best_pair) == [["A"], ["B"]]
    for name in ("A", "B", "C"):
        stat = rec[name].statistic
        assert stat is not None and stat > 0.9
        assert not rec[name].flagged
    assert set(result.statistics()) == {"A", "B", "C", "D"}


def test_same_seed_gives_the_same_result(data: AnnData, result: DoubletPopulationResult) -> None:
    again = doublet_population_test(data, data.obs["truth"].astype(str), seed=0, auroc_max=0.80)
    assert again == result


def test_small_clusters_and_too_few_groups_are_not_tested(data: AnnData) -> None:
    truth = data.obs["truth"].astype(str)
    keep = truth.isin(["A", "D"]).to_numpy()
    res = doublet_population_test(data[keep].copy(), truth[keep], seed=0, auroc_max=0.80)
    assert [(r.cluster, r.statistic, r.flagged) for r in res.records] == [
        ("A", None, False),
        ("D", None, False),
    ]
    res = doublet_population_test(data, truth, seed=0, auroc_max=0.80, min_cluster_size=200)
    assert res.records == [] and res.flagged == []


def _run(data: AnnData, cfg: GridConfig, params: PipelineParams):
    return Pipeline(data.copy(), cfg).run(params)


def test_pipeline_removes_the_doublet_cluster(data: AnnData) -> None:
    rr = _run(data, GridConfig(), PARAMS)
    truth = data.obs["truth"]
    assert (truth.loc[rr.clusters.index] != "D").all()
    assert (truth == "A").sum() == (truth.loc[rr.clusters.index] == "A").sum()
    info = rr.effective["doublet_clusters"]
    assert isinstance(info, dict)
    assert len(info["removed"]) == 1 and info["n_cells"] > 0
    assert info["statistics"][info["removed"][0]] < 0.80
    parents = info["parents"][info["removed"][0]]
    assert len(parents) == 2 and all(parents)
    assert rr.n_input == rr.n_cells + rr.n_removed_qc + rr.n_removed_doublet
    assert rr.doublet_scores.index.equals(rr.clusters.index)
    assert rr.n_clusters == rr.clusters.nunique() == 3
    assert any("doublet populations" in w for w in rr.warnings)
    again = _run(data, GridConfig(), PARAMS)
    assert again.clusters.equals(rr.clusters)
    assert again.effective["doublet_clusters"] == info


def test_pipeline_keeps_the_block_without_doublet_removal(data: AnnData) -> None:
    rr = _run(data, GridConfig(), PARAMS.model_copy(update={"doublet_removal": False}))
    assert (data.obs["truth"].loc[rr.clusters.index] == "D").sum() == N_DOUBLET
    assert "doublet_clusters" not in rr.effective
    assert rr.n_removed_doublet == 0


def test_config_switch_disables_the_step(data: AnnData) -> None:
    rr = _run(data, GridConfig(doublet_cluster_removal=False), PARAMS)
    assert "doublet_clusters" not in rr.effective
    assert (data.obs["truth"].loc[rr.clusters.index] == "D").sum() > N_DOUBLET // 2
    assert not any("doublet populations" in w for w in rr.warnings)


def test_config_defaults() -> None:
    cfg = GridConfig()
    assert cfg.doublet_cluster_removal is True
    assert cfg.doublet_cluster_auroc_max == 0.80
