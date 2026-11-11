"""Audit mode on a processed object, and the CLI entry points."""

from __future__ import annotations

import warnings

import pytest
import scanpy as sc
from typer.testing import CliRunner

from recluster import __version__
from recluster.cli import app
from recluster.data import DataError
from recluster.engine import run_claim
from recluster.spec import ClaimSpec


@pytest.fixture(scope="module")
def processed_path(sim, tmp_path_factory):
    a = sim.copy()
    a.layers["counts"] = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    a.obs["tool_label"] = a.obs["cell_type"].astype(str)
    del a.obs["cell_type"]
    p = tmp_path_factory.mktemp("proc") / "processed.h5ad"
    a.write_h5ad(p)
    return p


def _label_only_spec(claim_type="population_exists", labels="tool_label"):
    d = {
        "claim_type": claim_type,
        "name": "Mono as labelled",
        "columns": {"sample": ["donor", "condition"], "batch": "donor", "labels": labels},
        "target": {"label_value": "Mono"},
    }
    if claim_type == "composition_shift":
        d["columns"]["condition"] = "condition"
        d["comparison"] = {"group_a": "stim", "group_b": "ctrl", "direction": "higher_in_a"}
    return ClaimSpec.model_validate(d)


def test_audit_uses_stored_counts_and_label_matching(processed_path, quick_grid, tmp_path):
    res = run_claim(
        processed_path,
        _label_only_spec(),
        out_dir=tmp_path,
        grid=quick_grid,
        audit=True,
        progress=None,
    )
    assert res.mode == "audit"
    assert res.data_summary["counts_source"] == "layers['counts']"
    assert res.runs[0].match["mode"] == "labels"
    assert res.verdict.label == "ROBUST"
    assert res.as_labelled["n_cells"] > 0


def test_audit_composition_reports_as_labelled(processed_path, quick_grid, tmp_path):
    res = run_claim(
        processed_path,
        _label_only_spec("composition_shift"),
        out_dir=tmp_path,
        grid=quick_grid,
        audit=True,
        progress=None,
    )
    assert res.verdict.label == "ROBUST"
    assert res.as_labelled["holds"] is True
    assert res.as_labelled["n_samples_a"] == 6


def test_audit_refuses_raw_counts_in_x(sim_path, quick_grid, tmp_path):
    with pytest.raises(DataError, match="audit mode expects"):
        run_claim(
            sim_path,
            _label_only_spec(labels="cell_type"),
            out_dir=tmp_path,
            grid=quick_grid,
            audit=True,
            progress=None,
        )


def test_cli_version_and_run(sim_path, tmp_path):
    runner = CliRunner()
    assert runner.invoke(app, ["version"]).output.strip() == __version__
    claim = tmp_path / "claim.yaml"
    claim.write_text(
        "claim_type: population_exists\nname: cli smoke\n"
        "columns: {sample: [donor, condition], labels: cell_type}\n"
        "target: {markers: [Mono_marker_0, Mono_marker_1, Mono_marker_2], label_value: Mono}\n"
    )
    r = runner.invoke(
        app,
        [
            "run",
            "--data",
            str(sim_path),
            "--claim",
            str(claim),
            "--out",
            str(tmp_path),
            "--max-runs",
            "2",
            "-q",
        ],
    )
    assert r.exit_code == 0, r.output
    assert r.output.startswith("ROBUST\t2/2 runs")
    assert (tmp_path / "cli-smoke" / "report.html").exists()


def test_cli_reports_bad_spec_cleanly(sim_path, tmp_path):
    claim = tmp_path / "bad.yaml"
    claim.write_text("claim_type: composition_shift\nname: x\ntarget: {markers: [A]}\n")
    r = CliRunner().invoke(app, ["run", "--data", str(sim_path), "--claim", str(claim)])
    assert r.exit_code == 2
    assert "invalid claim spec" in r.output
