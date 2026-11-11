"""Shared fixtures: one synthetic dataset per session, written once to disk."""

from __future__ import annotations

from pathlib import Path

import pytest
from anndata import AnnData

from recluster.config import GridAxes, GridConfig, Thresholds
from recluster.simulate import SimConfig, simulate


@pytest.fixture(scope="session")
def sim() -> AnnData:
    return simulate(SimConfig(seed=0))


@pytest.fixture(scope="session")
def sim_path(sim: AnnData, tmp_path_factory: pytest.TempPathFactory) -> Path:
    p = tmp_path_factory.mktemp("data") / "synthetic.h5ad"
    sim.write_h5ad(p)
    return p


@pytest.fixture(scope="session")
def truth(sim: AnnData) -> dict:
    return dict(sim.uns["truth"])


@pytest.fixture
def quick_grid() -> GridConfig:
    """Baseline + resolution 0.5 + seed 1 + one random combination = 4 runs."""
    return GridConfig(
        axes=GridAxes(
            min_genes=[200],
            max_mito_pct=[10.0],
            doublet_removal=[True],
            normalization=["log_total"],
            n_hvgs=[2000],
            n_pcs=[30],
            resolution=[0.5, 1.0],
            batch_correction=[True],
            seed=[0, 1],
        ),
        n_random=1,
    )


@pytest.fixture
def thresholds() -> Thresholds:
    return Thresholds()


def mono_markers(truth: dict) -> list[str]:
    return list(truth["markers"]["Mono"][:5])
