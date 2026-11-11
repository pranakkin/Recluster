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
