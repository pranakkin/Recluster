"""Grid construction and verdict rules."""

from __future__ import annotations

from recluster.config import GridConfig, PipelineParams, Thresholds
from recluster.grid import build_grid, describe_grid
from recluster.results import RunEvaluation
from recluster.verdict import FRAGILE, NOT_SUPPORTED, ROBUST, compute_verdict


def test_grid_structure_default():
    cfg = GridConfig()
    runs = build_grid(cfg, has_batch=True)
    assert runs[0] == cfg.baseline
    assert len(runs) == len(set(runs))
    d = describe_grid(runs)
    assert d["one_at_a_time"] == 17 and d["random"] == 6 and d["total"] == 24
    for r in runs[1 : 1 + d["one_at_a_time"]]:
        assert len(r.diff(runs[0])) == 1


def test_grid_drops_batch_axis_without_batch_column():
    runs = build_grid(GridConfig(), has_batch=False)
    assert all(r.batch_correction is False for r in runs)
    assert describe_grid(runs)["one_at_a_time"] == 16


def test_grid_is_deterministic_and_capped():
    a = build_grid(GridConfig(random_seed=3), has_batch=True)
    b = build_grid(GridConfig(random_seed=3), has_batch=True)
    assert a == b
    assert build_grid(GridConfig(random_seed=4), has_batch=True) != a
    assert len(build_grid(GridConfig(max_runs=5), has_batch=True)) == 5


def _ev(i, params, holds, cell=None, base=None):
    base = base or PipelineParams()
    return RunEvaluation(
        index=i,
        label=params.label(base),
        params=params.model_dump(),
        is_baseline=(i == 0),
        is_one_at_a_time=(len(params.diff(base)) == 1),
        n_cells=1000,
        n_removed_qc=0,
        n_removed_doublet=0,
        n_clusters=5,
        runtime_s=1.0,
        warnings=[],
        match={},
        holds=holds,
        holds_cell_level=cell,
        replicate_ok=holds,
        recovered=holds,
        key_stats={},
    )


def test_verdict_robust_fragile_not_supported():
    runs = build_grid(GridConfig(n_random=0), has_batch=False)
    base = runs[0]
    thr = Thresholds()
    allgood = [_ev(i, p, True, base=base) for i, p in enumerate(runs)]
    assert compute_verdict(allgood, base, thr, "population_exists").label == ROBUST
    half = [_ev(i, p, i % 2 == 0, base=base) for i, p in enumerate(runs)]
    v = compute_verdict(half, base, thr, "population_exists")
    assert v.label == FRAGILE and v.breakers and v.oat_flips
    none = [_ev(i, p, False, base=base) for i, p in enumerate(runs)]
    assert compute_verdict(none, base, thr, "population_exists").label == NOT_SUPPORTED


def test_verdict_attributes_failure_to_the_changed_axis():
    runs = build_grid(GridConfig(n_random=0), has_batch=False)
    base = runs[0]
    evs = [_ev(i, p, p.resolution != 0.3, base=base) for i, p in enumerate(runs)]
    v = compute_verdict(evs, base, Thresholds(), "population_exists")
    assert v.baseline_holds
    assert v.oat_flips == [
        {"label": "resolution=0.3", "holds": False, "params_changed": {"resolution": 0.3}}
    ]
    top = v.breakers[0]
    assert top["axis"] == "resolution" and top["value"] == 0.3 and top["n_holds"] == 0


def test_pseudoreplication_flag():
    runs = build_grid(GridConfig(n_random=0), has_batch=False)
    base = runs[0]
    evs = [_ev(i, p, False, cell=True, base=base) for i, p in enumerate(runs)]
    v = compute_verdict(evs, base, Thresholds(), "differential_expression")
    assert v.pseudoreplication_flag and v.label == NOT_SUPPORTED
    assert "pseudoreplication" in v.summary
    v2 = compute_verdict(evs, base, Thresholds(), "population_exists")
    assert not v2.pseudoreplication_flag
