"""End-to-end through the engine on a 4-run grid: outputs, verdicts, determinism."""

from __future__ import annotations

import csv
import json

import pytest

from recluster.engine import run_claim
from recluster.report.html import render_report
from recluster.spec import ClaimSpec

from .conftest import mono_markers

COLS = {
    "condition": "condition",
    "sample": ["donor", "condition"],
    "batch": "donor",
    "labels": "cell_type",
}


def _population_spec(truth):
    return ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": "Mono exists",
            "columns": COLS,
            "target": {"markers": mono_markers(truth), "label_value": "Mono"},
        }
    )


def _de_spec(truth, genes, direction, name):
    return ClaimSpec.model_validate(
        {
            "claim_type": "differential_expression",
            "name": name,
            "columns": COLS,
            "target": {"markers": mono_markers(truth), "label_value": "Mono"},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": direction},
            "genes": genes,
        }
    )


@pytest.fixture(scope="module")
def population_result(sim_path, truth, tmp_path_factory):
    from recluster.config import GridAxes, GridConfig

    grid = GridConfig(
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
    out = tmp_path_factory.mktemp("out")
    res = run_claim(sim_path, _population_spec(truth), out_dir=out, grid=grid, progress=None)
    return res, out


def test_population_claim_is_robust_and_writes_outputs(population_result):
    res, out = population_result
    assert res.verdict.label == "ROBUST"
    assert len(res.runs) == 4 and res.runs[0].is_baseline
    folder = out / "mono-exists"
    assert (folder / "report.html").exists()
    prov = json.loads((folder / "provenance.json").read_text())
    assert prov["verdict"]["label"] == "ROBUST"
    assert prov["input"]["sha256"] and prov["environment"]["packages"]["scanpy"]
    assert len(prov["runs"]) == 4
    with open(folder / "runs.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 4 and rows[0]["label"] == "baseline"
    assert "labelled `Mono`" in res.reference_source


def test_report_is_self_contained(population_result):
    res, _ = population_result
    html = render_report(res)
    assert "ROBUST" in html and "<svg" in html and "Methods" in html
    # no external assets: no scripts, stylesheets, fonts or images fetched from anywhere
    assert "<script" not in html and "<link" not in html and "@import" not in html
    assert 'src="http' not in html and "url(http" not in html


def test_false_de_claim_not_supported(sim_path, truth, quick_grid, tmp_path):
    spec = _de_spec(truth, list(truth["null_genes"][:2]), "either", "nulls change")
    res = run_claim(sim_path, spec, out_dir=tmp_path, grid=quick_grid, progress=None)
    assert res.verdict.label == "NOT SUPPORTED"
    assert res.as_labelled is not None and res.as_labelled["holds"] is False


def test_true_de_claim_robust(sim_path, truth, quick_grid, tmp_path):
    spec = _de_spec(truth, list(truth["de_genes"][:2]), "up_in_a", "planted up")
    res = run_claim(sim_path, spec, out_dir=tmp_path, grid=quick_grid, progress=None)
    assert res.verdict.label == "ROBUST"
    assert res.as_labelled["holds"] is True


def test_determinism_same_config_same_numbers(sim_path, truth, quick_grid, tmp_path):
    spec = _population_spec(truth)
    a = run_claim(sim_path, spec, out_dir=tmp_path / "a", grid=quick_grid, progress=None)
    b = run_claim(sim_path, spec, out_dir=tmp_path / "b", grid=quick_grid, progress=None)

    def strip(prov):
        runs = []
        for r in prov["runs"]:
            r = dict(r)
            r.pop("runtime_s")
            runs.append(r)
        return runs, prov["verdict"]

    assert strip(a.provenance) == strip(b.provenance)
    pa = json.loads((tmp_path / "a" / "mono-exists" / "provenance.json").read_text())
    pb = json.loads((tmp_path / "b" / "mono-exists" / "provenance.json").read_text())
    assert strip(pa) == strip(pb)


# ------------------------------------------------- determinism of comparison claims
# A differential-expression claim exercises the replicate-level DE test, the time-of-day check
# (Mode A with a collection-time column, Mode B without) and the as-labelled evaluation, none
# of which the population claim above touches.

EARLY_DONORS = ("donor0", "donor1", "donor2")
TIME_LABELS = ("NO TIME CONFOUND DETECTED", "POSSIBLE TIME CONFOUND", "NOT ASSESSABLE")


def _timed_sim_path(sim, tmp_path_factory):
    """The sim with one collection time per donor, written to a temporary h5ad.

    Donors 0-2 were collected at 09:00 and donors 3-5 at 15:30; both conditions of a donor
    share the time, so time is perfectly balanced between stim and ctrl.
    """
    a = sim.copy()
    early = a.obs["donor"].astype(str).isin(EARLY_DONORS)
    a.obs["collected"] = early.map({True: "09:00", False: "15:30"})
    path = tmp_path_factory.mktemp("timed") / "synthetic_timed.h5ad"
    a.write_h5ad(path)
    return path


@pytest.fixture(scope="module")
def timed_sim_path(sim, tmp_path_factory):
    return _timed_sim_path(sim, tmp_path_factory)


def _two_run_grid():
    """Baseline plus the first one-at-a-time run."""
    from recluster.config import GridConfig

    return GridConfig(max_runs=2, n_random=0)


def _timed_de_spec(truth):
    return ClaimSpec.model_validate(
        {
            "claim_type": "differential_expression",
            "name": "planted up timed",
            "columns": {**COLS, "collection_time": "collected"},
            "target": {"markers": mono_markers(truth), "label_value": "Mono"},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "up_in_a"},
            "genes": list(truth["de_genes"][:2]),
        }
    )


def _without_timing(prov):
    """The provenance minus the wall-clock fields, for a deep comparison of everything else."""
    out = {
        k: v for k, v in prov.items() if k not in ("started_at", "finished_at", "total_runtime_s")
    }
    out["runs"] = [{k: v for k, v in r.items() if k != "runtime_s"} for r in prov["runs"]]
    return out


def test_determinism_de_claim_with_times(timed_sim_path, truth, tmp_path):
    spec = _timed_de_spec(truth)
    grid = _two_run_grid()
    a = run_claim(timed_sim_path, spec, out_dir=tmp_path / "a", grid=grid, progress=None)
    b = run_claim(timed_sim_path, spec, out_dir=tmp_path / "b", grid=grid, progress=None)
    assert len(a.runs) == 2
    assert _without_timing(a.provenance) == _without_timing(b.provenance)

    prov = a.provenance
    tc = prov["time_check"]
    assert tc["mode"] == "A" and tc["label"] in TIME_LABELS
    # one time per donor and both conditions within every donor: time cannot track condition
    assert tc["label"] == "NO TIME CONFOUND DETECTED"
    assert tc["balance"]["imbalanced"] is False
    assert tc["balance"]["n_units_a"] == 6 and tc["balance"]["n_units_b"] == 6
    assert tc["adjustment"] is None  # only an imbalanced design is refitted
    assert prov["verdict"]["time_check_label"] == "NO TIME CONFOUND DETECTED"
    assert prov["verdict"]["label_before_time_check"] == prov["verdict"]["label"] == "ROBUST"
    assert prov["as_labelled"]["holds"] is True
    assert isinstance(prov["merged_runs"], list)

    on_disk = json.loads((tmp_path / "a" / "planted-up-timed" / "provenance.json").read_text())
    assert _without_timing(on_disk) == _without_timing(prov)


def test_determinism_de_claim_without_times(sim_path, truth, tmp_path):
    spec = _de_spec(truth, list(truth["de_genes"][:2]), "up_in_a", "planted up")
    grid = _two_run_grid()
    a = run_claim(sim_path, spec, out_dir=tmp_path / "a", grid=grid, progress=None)
    b = run_claim(sim_path, spec, out_dir=tmp_path / "b", grid=grid, progress=None)
    assert len(a.runs) == 2
    assert _without_timing(a.provenance) == _without_timing(b.provenance)

    prov = a.provenance
    tc = prov["time_check"]
    assert tc["mode"] == "B" and tc["label"] == "NOT ASSESSABLE"  # the sim has no clock genes
    assert tc["panel"]["assessable"] is False and tc["panel"]["genes_usable"] == []
    assert tc["balance"] is None and tc["adjustment"] is None
    assert prov["verdict"]["time_check_label"] == "NOT ASSESSABLE"
    assert prov["verdict"]["label_before_time_check"] == prov["verdict"]["label"] == "ROBUST"
    assert prov["as_labelled"]["holds"] is True
    assert isinstance(prov["merged_runs"], list)

    on_disk = json.loads((tmp_path / "a" / "planted-up" / "provenance.json").read_text())
    assert _without_timing(on_disk) == _without_timing(prov)
