"""First-use CLI behaviour: `recluster check`, readable file errors, help and version."""

from __future__ import annotations

import re
from dataclasses import dataclass

import pytest
from typer.testing import CliRunner

from recluster import __version__
from recluster.cli import app

OK_CLAIM = (
    "claim_type: population_exists\nname: check smoke\n"
    "columns: {sample: [donor, condition], labels: cell_type}\n"
    "target: {markers: [Mono_marker_0, Mono_marker_1], label_value: Mono}\n"
)


@dataclass
class _Result:
    exit_code: int
    output: str


def _invoke(*args: object) -> _Result:
    r = CliRunner().invoke(app, [str(a) for a in args], env={"COLUMNS": "200"})
    # CI terminals get colour codes in rich help output
    return _Result(r.exit_code, re.sub(r"\x1b\[[0-9;]*m", "", r.output))


def _claim(tmp_path, text: str, name: str = "claim.yaml"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_check_accepts_a_valid_claim(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    r = _invoke("check", "--data", sim_path, "--claim", claim)
    assert r.exit_code == 0, r.output
    assert "OK: 'check smoke' matches" in r.output
    assert "raw counts in .X" in r.output
    assert f"Next: recluster run --data {sim_path} --claim {claim}" in r.output


def test_check_reports_a_misspelt_column(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM.replace("labels: cell_type", "labels: cel_type"))
    r = _invoke("check", "--data", sim_path, "--claim", claim)
    assert r.exit_code == 2
    assert "not in .obs" in r.output
    assert "Did you mean 'cell_type'" in r.output


def test_check_audit_requires_labels(sim_path, tmp_path):
    claim = _claim(
        tmp_path,
        "claim_type: population_exists\nname: x\ntarget: {markers: [Mono_marker_0]}\n",
    )
    r = _invoke("check", "--data", sim_path, "--claim", claim, "--audit")
    assert r.exit_code == 2
    assert "audit mode requires" in r.output


@pytest.mark.parametrize(
    "text, expected",
    [
        (
            OK_CLAIM.replace("population_exists", "population_exist"),
            "Did you mean 'population_exists'",
        ),
        (OK_CLAIM.replace("markers:", "marker:"), "Did you mean 'markers'"),
        ('claim_type: population_exists\nname: "x\n', "not valid YAML at line"),
        ("", "is empty"),
        ("claim_type: population_exists\ntarget: {markers: [A]}\n", "`name` is required"),
    ],
)
def test_bad_claim_files_are_explained(sim_path, tmp_path, text, expected):
    claim = _claim(tmp_path, text)
    r = _invoke("check", "--data", sim_path, "--claim", claim)
    assert r.exit_code == 2
    assert r.output.startswith(f"invalid claim spec {claim}:")
    assert expected in r.output
    assert "Examples to copy: examples/claims/" in r.output


def test_swapped_data_and_claim(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    r = _invoke("check", "--data", claim, "--claim", sim_path)
    assert r.exit_code == 2
    assert "swapped" in r.output


def test_yaml_given_as_data(tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    r = _invoke("check", "--data", claim, "--claim", claim)
    assert r.exit_code == 2
    assert "could not be read as an AnnData" in r.output


def test_bad_grid_file_is_explained(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    grid = _claim(tmp_path, "axes: {resolutions: [0.5, 1.0]}\n", "grid.yaml")
    r = _invoke("run", "--data", sim_path, "--claim", claim, "--grid", grid)
    assert r.exit_code == 2
    assert r.output.startswith(f"invalid grid file {grid}:")
    assert "Did you mean 'resolution'" in r.output


def test_bad_thresholds_file_is_explained(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    thr = _claim(tmp_path, "robust_min: 0.9\n", "thr.yaml")
    r = _invoke("run", "--data", sim_path, "--claim", claim, "--thresholds", thr)
    assert r.exit_code == 2
    assert "invalid thresholds file" in r.output
    assert "Did you mean 'robust_min_survival'" in r.output


def test_run_prints_the_report_path(sim_path, tmp_path):
    claim = _claim(tmp_path, OK_CLAIM)
    out = tmp_path / "out"
    r = _invoke("run", "--data", sim_path, "--claim", claim, "--out", out, "--max-runs", "2", "-q")
    assert r.exit_code == 0, r.output
    first, second = r.output.splitlines()[:2]
    assert first.startswith("ROBUST\t2/2 runs")
    report = out / "check-smoke" / "report.html"
    assert second == f"report: {report}"
    assert report.exists()


def test_version_flag_and_help():
    assert _invoke("--version").output.strip() == __version__
    assert _invoke("version").output.strip() == __version__
    r = _invoke("--help")
    assert r.exit_code == 0
    assert "check" in r.output
    assert "ROBUST" in r.output
    assert "Start here" in r.output


def test_bench_rejects_unknown_part():
    r = _invoke("bench", "--part", "nope")
    assert r.exit_code == 2


def test_demo_help_lists_its_options():
    r = _invoke("demo", "--help")
    assert r.exit_code == 0
    assert "--no-open" in r.output
    assert "--jobs" in r.output
