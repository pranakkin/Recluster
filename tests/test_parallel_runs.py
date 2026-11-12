"""Grid runs in worker processes give the same results as running them one after another."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from recluster.config import GridConfig
from recluster.engine import run_claim, slugify
from recluster.spec import ClaimSpec

from .conftest import mono_markers

COLS = {
    "condition": "condition",
    "sample": ["donor", "condition"],
    "batch": "donor",
    "labels": "cell_type",
}


def _de_spec(truth) -> ClaimSpec:
    return ClaimSpec.model_validate(
        {
            "claim_type": "differential_expression",
            "name": "parallel de",
            "columns": COLS,
            "target": {"markers": mono_markers(truth), "label_value": "Mono"},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "up_in_a"},
            "genes": list(truth["de_genes"][:2]),
        }
    )


def _marker_population_spec(truth) -> ClaimSpec:
    # no label_value: the reference comes from the first run whose marker match qualifies
    return ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": "parallel population",
            "columns": {k: v for k, v in COLS.items() if k != "labels"},
            "target": {"markers": mono_markers(truth)},
        }
    )


def _composition_spec(truth) -> ClaimSpec:
    return ClaimSpec.model_validate(
        {
            "claim_type": "composition_shift",
            "name": "parallel composition",
            "columns": COLS,
            "target": {"markers": mono_markers(truth), "label_value": "Mono"},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "higher_in_a"},
        }
    )


SPECS = {
    "de": _de_spec,
    "marker_population": _marker_population_spec,
    "composition": _composition_spec,
}


def _comparable(prov: dict) -> dict:
    """Provenance without wall-clock fields and the execution block."""
    skip = ("started_at", "finished_at", "total_runtime_s", "execution")
    out = {k: v for k, v in prov.items() if k not in skip}
    out["runs"] = [{k: v for k, v in r.items() if k != "runtime_s"} for r in prov["runs"]]
    return out


def _runs_csv(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r.pop("runtime_s")
    return rows


@pytest.fixture(scope="module", params=sorted(SPECS))
def both(request, sim_path, truth, tmp_path_factory):
    spec = SPECS[request.param](truth)
    grid = GridConfig(max_runs=4, n_random=0)
    out = {}
    for jobs in (1, 2):
        folder = tmp_path_factory.mktemp(f"{request.param}_jobs{jobs}")
        res = run_claim(sim_path, spec, out_dir=folder, grid=grid, progress=None, jobs=jobs)
        out[jobs] = (res, folder / slugify(spec.name))
    return out


def test_same_provenance(both):
    (seq, seq_dir), (par, par_dir) = both[1], both[2]
    assert len(seq.runs) == 4
    assert _comparable(par.provenance) == _comparable(seq.provenance)
    on_disk_seq = json.loads((seq_dir / "provenance.json").read_text())
    on_disk_par = json.loads((par_dir / "provenance.json").read_text())
    assert _comparable(on_disk_par) == _comparable(on_disk_seq)
    assert seq.verdict.to_dict() == par.verdict.to_dict()
    assert seq.reference_source == par.reference_source


def test_same_runs_csv(both):
    (_, seq_dir), (_, par_dir) = both[1], both[2]
    assert _runs_csv(par_dir / "runs.csv") == _runs_csv(seq_dir / "runs.csv")


def test_execution_recorded(both):
    assert both[1][0].provenance["execution"] == {"jobs": 1, "workers": 1}
    assert both[2][0].provenance["execution"] == {"jobs": 2, "workers": 2}


def test_jobs_must_be_positive(sim_path, truth):
    with pytest.raises(ValueError, match="jobs"):
        run_claim(sim_path, _de_spec(truth), progress=None, write=False, jobs=0)
