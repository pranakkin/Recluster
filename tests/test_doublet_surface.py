"""Where the doublet-population step shows up: methods text, report, runs.csv, provenance, and
the planted doublet cases of the benchmark (built here on synthetic data, never run)."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from recluster.bench import runner
from recluster.bench.planted import plant_doublets, plant_held_out_doublets
from recluster.config import GridAxes, GridConfig, PipelineParams, Thresholds
from recluster.engine import run_claim
from recluster.report.html import _doublet_population_html, render_report
from recluster.report.methods import methods_paragraph
from recluster.results import RunEvaluation
from recluster.spec import ClaimSpec

# ------------------------------------------------------------------- methods text


def _methods(grid: GridConfig) -> str:
    spec = ClaimSpec.model_validate(
        {"claim_type": "population_exists", "name": "m", "target": {"markers": ["G1"]}}
    )
    data = {"n_cells": 100, "n_genes": 50, "counts_source": ".X", "sha256": "0" * 64}
    counts = {"one_at_a_time": 1, "random": 0}
    return methods_paragraph(spec, grid, Thresholds(), data, 2, counts, "none")


def test_methods_describe_the_doublet_population_step() -> None:
    text = _methods(GridConfig())
    assert "every cluster of at least 10 cells was compared with simulated doublets" in text
    assert "cross-validated logistic-regression classifier" in text
    assert "below 0.80 for some pair" in text
    assert "removed as a doublet population" in text
    assert "below 0.75" in _methods(GridConfig(doublet_cluster_auroc_max=0.75))
    off = _methods(GridConfig(doublet_cluster_removal=False))
    assert "not tested as doublet populations (doublet_cluster_removal off)" in off
    assert "cross-validated" not in off


# ------------------------------------------------- engine run with a doublet block

N_REAL, N_DOUBLET = 150, 80


@pytest.fixture(scope="module")
def doublet_data_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Three Poisson populations and a block of A+B doublets (as in test_doublet_population)."""
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
    path = tmp_path_factory.mktemp("dbl") / "doublet_block.h5ad"
    AnnData(X=sp.csr_matrix(x), obs=obs, var=var).write_h5ad(path)
    return path


@pytest.fixture(scope="module")
def engine_out(doublet_data_path: Path, tmp_path_factory: pytest.TempPathFactory):
    """Baseline (doublet removal on) and the one-at-a-time run with it off."""
    grid = GridConfig(
        baseline=PipelineParams(
            min_genes=10, resolution=0.5, n_pcs=10, batch_correction=False, doublet_removal=True
        ),
        axes=GridAxes(
            min_genes=[10],
            max_mito_pct=[10.0],
            doublet_removal=[False, True],
            normalization=["log_total"],
            n_hvgs=[2000],
            n_pcs=[10],
            resolution=[0.5],
            batch_correction=[False],
            seed=[0],
        ),
        n_random=0,
    )
    spec = ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": "A exists",
            "columns": {"labels": "truth"},
            "target": {"label_value": "A", "min_cells": 20},
        }
    )
    out = tmp_path_factory.mktemp("out")
    res = run_claim(doublet_data_path, spec, out_dir=out, grid=grid, progress=None)
    return res, out / "a-exists"


def test_engine_records_the_removed_cluster(engine_out) -> None:
    res, folder = engine_out
    assert len(res.runs) == 2
    base, off = res.runs
    assert base.params["doublet_removal"] and not off.params["doublet_removal"]
    info = base.effective_params["doublet_clusters"]
    assert len(info["removed"]) == 1 and info["n_cells"] > 0
    removed = info["removed"][0]
    assert info["statistics"][removed] < 0.80
    assert "doublet_clusters" not in off.effective_params
    # the matched population (A) carries its own statistic, well above the threshold
    pop = base.match["doublet_population"]
    assert pop["auroc_max"] == 0.80
    assert pop["statistics"] and all(v > 0.9 for v in pop["statistics"].values())
    assert set(pop["statistics"]) == set(base.match["clusters"])
    assert off.match["doublet_population"] is None

    with open(folder / "runs.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert int(rows[0]["n_removed_doublet_cluster"]) == info["n_cells"]
    assert int(rows[0]["n_removed_doublet"]) >= info["n_cells"]
    assert rows[0]["doublet_clusters_removed"].startswith(f"{removed}:0.")
    assert rows[1]["n_removed_doublet_cluster"] == "" and rows[1]["doublet_clusters_removed"] == ""

    prov = json.loads((folder / "provenance.json").read_text())
    assert prov["runs"][0]["effective_params"]["doublet_clusters"]["removed"] == [removed]
    assert prov["runs"][0]["n_removed_doublet"] == base.n_removed_doublet
    assert prov["grid"]["doublet_cluster_auroc_max"] == 0.80


def test_report_shows_the_removed_cluster_block(engine_out) -> None:
    res, folder = engine_out
    html = render_report(res)
    removed = res.runs[0].effective_params["doublet_clusters"]["removed"][0]
    assert "Doublet-population removal (baseline run)" in html
    assert "Matched population, AUROC vs simulated doublets" in html
    assert "<th>Removed cluster</th>" in html
    assert f"<td class='mono'>{removed}</td>" in html
    parents = res.runs[0].effective_params["doublet_clusters"]["parents"][removed]
    assert " + ".join(",".join(side) for side in parents) in html
    assert "in 1 of the 1 runs where the step ran (run 0: cluster" in html
    assert "removed as a doublet population" in html  # methods paragraph
    assert 'doublet_cluster_auroc_max</span> = <span class="mono">0.8</span>' in html
    assert "Doublet-population removal (baseline run)" in (folder / "report.html").read_text()


def _run(effective: dict[str, Any], match: dict[str, Any]) -> RunEvaluation:
    return RunEvaluation(
        index=0,
        label="baseline",
        params={"doublet_removal": True},
        is_baseline=True,
        is_one_at_a_time=False,
        n_cells=100,
        n_removed_qc=0,
        n_removed_doublet=30,
        n_clusters=3,
        runtime_s=0.0,
        warnings=[],
        match=match,
        holds=True,
        holds_cell_level=None,
        replicate_ok=None,
        recovered=None,
        key_stats={},
        effective_params=effective,
    )


def test_block_shows_parent_groups_and_untested_clusters() -> None:
    info = {
        "removed": ["7"],
        "n_cells": 25,
        "auroc_max": 0.8,
        "statistics": {"2": None, "7": 0.61},
        "parents": {"7": [["1", "4"], ["0"]]},
    }
    match = {"clusters": ["2", "9"], "doublet_population": {"statistics": {"2": None}}}
    base = _run({"doublet_clusters": info}, match)
    html = _doublet_population_html(SimpleNamespace(runs=[base]), base)  # type: ignore[arg-type]
    assert "1,4 + 0" in html and "0.610" in html and "1 (25)" in html
    assert "cluster 2: not computed (fewer than two other cell-type groups)" in html
    assert "cluster 9: not tested (under 10 cells)" in html
    off = _run({}, {"clusters": []})
    html = _doublet_population_html(SimpleNamespace(runs=[off]), off)  # type: ignore[arg-type]
    assert "Not run in the baseline run" in html


# ------------------------------------------------------- planted doublet cases

TYPES = ("NK", "CD14 Mono", "CD8 T", "B")
N_PER_TYPE = 160


@pytest.fixture(scope="module")
def tagged() -> tuple[AnnData, pd.Series]:
    """One gene per cell, so every row of a doublet names its two source cells."""
    n = N_PER_TYPE * len(TYPES)
    x = sp.identity(n, dtype=np.float32, format="csr") * 3
    names = [f"cell{i}" for i in range(n)]
    adata = AnnData(X=x, var=pd.DataFrame(index=[f"g{i}" for i in range(n)]))
    adata.obs_names = names
    ref = pd.Series(np.repeat(TYPES, N_PER_TYPE), index=names, name="reference")
    return adata, ref


def _sources(case_adata: AnnData, n_orig: int) -> np.ndarray:
    """(n_doublets, 2) source-cell indices of the appended rows."""
    x = sp.csr_matrix(case_adata.X[n_orig:])
    assert (np.diff(x.indptr) == 2).all() and (x.data == 3).all()
    return x.indices.reshape(-1, 2)


def test_plant_doublets_uses_the_requested_pair_and_count(tagged) -> None:
    adata, ref = tagged
    (case,) = plant_doublets(adata, ref, type_a="NK", type_b="CD14 Mono", n=40, seed=11)
    out = case.adata
    assert out.n_obs == adata.n_obs + 40
    assert int(out.obs["injected_doublet"].sum()) == 40
    assert not out.obs["injected_doublet"].iloc[: adata.n_obs].any()
    src = np.sort(_sources(out, adata.n_obs), axis=1)
    types = ref.to_numpy()
    assert (types[src[:, 0]] == "NK").all() and (types[src[:, 1]] == "CD14 Mono").all()
    assert len(set(src[:, 0])) == 40 and len(set(src[:, 1])) == 40  # no cell used twice
    again = plant_doublets(adata, ref, type_a="NK", type_b="CD14 Mono", n=40, seed=11)[0]
    assert (again.adata.X != out.X).nnz == 0
    with pytest.raises(ValueError, match="'NK' cells"):
        plant_doublets(adata, ref, type_a="NK", type_b="B", n=N_PER_TYPE + 1)


def test_default_doublet_case_is_unchanged(tagged) -> None:
    adata, ref = tagged
    (case,) = plant_doublets(adata, ref, seed=0)
    assert case.key == "doublets" and case.extra == {"n_injected": 150}
    assert case.spec.name == (
        "Planted doublets: a CD14+ LYZ+ MS4A1+ CD79A+ population exists (PBMC3k)"
    )
    assert case.spec.target.markers == ["CD14", "LYZ", "MS4A1", "CD79A"]
    assert case.title == "Injected synthetic doublets forming a fake population"
    assert "a random CD14 monocyte and a random B cell (seed 0)" in case.how_built
    assert "both monocyte and B-cell markers" in case.how_built
    src = np.sort(_sources(case.adata, adata.n_obs), axis=1)
    assert set(ref.to_numpy()[src.ravel()]) == {"CD14 Mono", "B"}


def _view(label: str, runs: list[tuple[bool, bool]]) -> SimpleNamespace:
    """A result with the given verdict and (holds, doublet_removal) per run."""
    return SimpleNamespace(
        verdict=SimpleNamespace(label=label, breakers=[]),
        runs=[SimpleNamespace(holds=h, params={"doublet_removal": d}) for h, d in runs],
    )


def test_held_out_cases_are_built_with_the_strict_rule(tagged) -> None:
    adata, ref = tagged
    dev = plant_doublets(adata, ref, seed=0)[0]
    h1, h2 = plant_held_out_doublets(adata, ref)
    expected = {
        "doublets_nk_mono": (
            ("NK", "CD14 Mono"),
            11,
            ["GNLY", "NKG7", "CD14", "LYZ"],
            "Planted doublets (held out): a GNLY+ NKG7+ CD14+ LYZ+ population exists (PBMC3k)",
        ),
        "doublets_cd8_b": (
            ("CD8 T", "B"),
            12,
            ["CD8A", "CD8B", "MS4A1", "CD79A"],
            "Planted doublets (held out): a CD8A+ CD8B+ MS4A1+ CD79A+ population exists (PBMC3k)",
        ),
    }
    for case in (h1, h2):
        pair, seed, markers, name = expected[case.key]
        assert case.spec.name == name
        assert case.spec.claim_type == "population_exists"
        assert case.spec.target.markers == markers and case.spec.target.min_cells == 20
        assert case.extra == {
            "n_injected": 150,
            "held_out": True,
            "parent_types": list(pair),
            "seed": seed,
        }
        assert f"(seed {seed})" in case.how_built and "Held-out case" in case.how_built
        assert case.caught_rule == dev.caught_rule and case.expected == dev.expected
        assert not case.audit and case.role == "planted"
        out = case.adata
        assert out.n_obs == adata.n_obs + 150
        assert int(out.obs["injected_doublet"].sum()) == 150
        src = np.sort(_sources(out, adata.n_obs), axis=1)
        got = ref.to_numpy()[src]
        assert {tuple(sorted(r)) for r in got} == {tuple(sorted(pair))}
        # strict rule: not ROBUST, held somewhere, and only with doublet removal off
        assert case.caught(_view("FRAGILE", [(True, False), (False, True)]))
        assert not case.caught(_view("ROBUST", [(True, False)]))
        assert not case.caught(_view("NOT SUPPORTED", [(False, False), (False, True)]))
        assert not case.caught(_view("FRAGILE", [(True, False), (True, True)]))


def test_runner_builds_the_held_out_jobs(tagged, tmp_path, monkeypatch) -> None:
    """Wiring only: the planted part builds jobs for the held-out cases; nothing is run."""
    adata, ref = tagged
    monkeypatch.setattr(runner, "load_pbmc3k", lambda _d: adata)
    monkeypatch.setattr(runner, "reference_labels", lambda _a: ref)
    monkeypatch.setattr(runner, "scrublet_recall_on_injected", lambda _a: {"recall": 0.0})
    for fn in ("plant_overclustering", "plant_qc_deletion", "plant_batch_effect"):
        monkeypatch.setattr(runner, fn, lambda *_a, **_k: [])
    jobs, meta, info = runner.build_planted_jobs(
        tmp_path, tmp_path / "out", GridConfig(), Thresholds()
    )
    assert [j["key"] for j in jobs] == ["doublets", "doublets_nk_mono", "doublets_cd8_b"]
    assert len({j["data_path"] for j in jobs}) == 3
    assert all(Path(j["data_path"]).exists() and not j["audit"] for j in jobs)
    assert meta["doublets_nk_mono"]["extra"]["held_out"] is True
    assert meta["doublets_cd8_b"]["rule"] == meta["doublets"]["rule"]
    assert set(info["scrublet_on_injected_held_out"]) == {"doublets_nk_mono", "doublets_cd8_b"}
