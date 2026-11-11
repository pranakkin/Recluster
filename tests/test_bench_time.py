"""Planted time-of-day cases: the construction, the claim specs, and the pre-registered rules.

As in test_bench.py, the real-data runs are not exercised. These tests cover what must be right
for the two time cases to mean anything: one unit per donor, rules that a tool finding nothing
cannot pass, and a TIME_SIGNAL provenance that matches the shipped day/night table.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from recluster.bench import runner
from recluster.bench.kang import (
    AFTERNOON_TIMES,
    MORNING_TIMES,
    TIME_CLAIM_GENE,
    TIME_SIGNAL,
    TIME_SIGNAL_S8_DISAGREE,
    KangCase,
    kang_cases,
    plant_time_signal,
    time_signal_day_night,
)
from recluster.bench.planted import plant_doublets
from recluster.bench.runner import _compact_time_check, _ResultView, results_table
from recluster.config import GridConfig, Thresholds
from recluster.spec import make_sample_key
from recluster.time_confound import load_rhythmic_table

DONORS = [f"patient_{i}" for i in range(8)]
MARKERS = ["CD14", "LYZ", "S100A8", "IL7R", "CD3D", "CCR7"]
FILLER = [f"filler_{i:02d}" for i in range(50)]


@pytest.fixture
def kang_small() -> AnnData:
    """Kang-shaped synthetic object: 8 donors x ctrl/stim x 2 cell types, Poisson counts."""
    rng = np.random.default_rng(0)
    genes = list(TIME_SIGNAL) + MARKERS + FILLER
    rows = []
    for d in DONORS:
        for label in ("ctrl", "stim"):
            for ct in ("CD14+ Monocytes", "CD4 T cells"):
                rows += [{"replicate": d, "label": label, "cell_type": ct}] * 10
    obs = pd.DataFrame(rows, index=[f"cell{i}" for i in range(len(rows))])
    for c in obs.columns:
        obs[c] = obs[c].astype("category")
    means = rng.gamma(2.0, 2.0, size=len(genes))
    x = rng.poisson(means[None, :], size=(len(obs), len(genes))).astype(np.float32)
    return AnnData(X=sp.csr_matrix(x), obs=obs, var=pd.DataFrame(index=genes))


def _dense(a: AnnData) -> np.ndarray:
    return np.asarray(a.X.todense()) if sp.issparse(a.X) else np.asarray(a.X)


def _time_cases(planted: AnnData) -> dict[str, KangCase]:
    cases = kang_cases(
        Path("kang.h5ad"),
        Path("split.h5ad"),
        Path("shuffle.h5ad"),
        {"A": DONORS[:4], "B": DONORS[4:]},
        Path("time.h5ad"),
        dict(planted.uns["planted_time"]["groups"]),
    )
    return {c.key: c for c in cases if c.role == "planted_time"}


def _summary(verdict: str, before: str | None, time_check: dict[str, Any] | None) -> dict[str, Any]:
    return {
        "verdict": verdict,
        "top_breakers": [],
        "runs": [{"holds": verdict != "NOT SUPPORTED", "params": {"doublet_removal": True}}],
        "label_before_time_check": before,
        "time_check": time_check,
    }


# ------------------------------------------------------------- construction


def test_plant_time_signal_is_deterministic_and_touches_only_clock_genes(kang_small):
    p0 = plant_time_signal(kang_small, seed=0)
    again = plant_time_signal(kang_small, seed=0)
    assert np.array_equal(_dense(p0), _dense(again))
    assert p0.uns["planted_time"] == again.uns["planted_time"]
    assert not np.array_equal(_dense(p0), _dense(plant_time_signal(kang_small, seed=1)))
    assert p0.n_obs == kang_small.n_obs and list(p0.var_names) == list(kang_small.var_names)
    before, after = _dense(kang_small), _dense(p0)
    changed = {
        g for j, g in enumerate(p0.var_names) if not np.array_equal(before[:, j], after[:, j])
    }
    assert changed == set(TIME_SIGNAL)
    assert "fake_group" not in kang_small.obs and "planted_time" not in kang_small.uns


def test_plant_time_signal_gives_one_time_and_one_group_per_donor(kang_small):
    p = plant_time_signal(kang_small, seed=0)
    per_donor = p.obs.groupby("replicate", observed=True)
    assert (per_donor["fake_group"].nunique() == 1).all()
    assert (per_donor["collection_time"].nunique() == 1).all()
    group = per_donor["fake_group"].first().astype(str)
    hours = per_donor["collection_time"].first().astype(float)
    assert group.value_counts().to_dict() == {"A": 4, "B": 4}
    assert sorted(hours[group == "A"]) == MORNING_TIMES
    assert sorted(hours[group == "B"]) == AFTERNOON_TIMES
    planted = p.uns["planted_time"]
    assert set(planted) == {
        "groups",
        "collection_time_by_donor",
        "factors_by_gene_and_donor",
        "genes",
    }
    assert planted["groups"] == group.to_dict()
    assert planted["collection_time_by_donor"] == hours.to_dict()
    assert planted["genes"] == sorted(TIME_SIGNAL)
    peak, amp = TIME_SIGNAL["PER1"]
    for d, t in planted["collection_time_by_donor"].items():
        expected = 2.0 ** (amp * math.cos(2 * math.pi * (t - peak) / 24.0))
        assert math.isclose(planted["factors_by_gene_and_donor"]["PER1"][d], expected)


# ----------------------------------------------------------------- the specs


def test_time_specs_use_one_sample_per_donor_and_validate(kang_small):
    planted = plant_time_signal(kang_small, seed=0)
    cases = _time_cases(planted)
    assert list(cases) == ["time_mode_b", "time_mode_a"]
    for c in cases.values():
        c.spec.check_against(planted)
        assert c.spec.columns.sample_key == ["replicate"]
        assert c.spec.columns.condition == "fake_group"
        assert c.spec.columns.labels == "cell_type"
        assert c.spec.genes == [TIME_CLAIM_GENE]
        assert c.spec.comparison is not None and c.spec.comparison.group_a == "B"  # afternoon
        assert "Each donor is one sample" in c.how_built and "8 units, not 16" in c.how_built
        assert make_sample_key(planted.obs, c.spec.columns.sample_key).nunique() == 8
    assert cases["time_mode_b"].spec.columns.collection_time is None
    assert cases["time_mode_a"].spec.columns.collection_time == "collection_time"


def test_time_cases_are_absent_without_a_time_path():
    cases = kang_cases(
        Path("kang.h5ad"), Path("split.h5ad"), Path("shuffle.h5ad"), {"A": ["p1"], "B": ["p2"]}
    )
    assert not [c for c in cases if c.role == "planted_time"]


# ----------------------------------------------------------------- the rules


def test_mode_b_rule_needs_the_naive_call_the_flag_and_a_downgrade(kang_small):
    rule = _time_cases(plant_time_signal(kang_small, seed=0))["time_mode_b"].passes
    possible = {"mode": "B", "label": "POSSIBLE TIME CONFOUND"}
    assert rule(_ResultView(_summary("FRAGILE", "ROBUST", possible)))
    assert rule(_ResultView(_summary("NOT SUPPORTED", "ROBUST", possible)))
    # the replicate-aware analysis did not call the planted effect: nothing to catch
    assert not rule(_ResultView(_summary("FRAGILE", "FRAGILE", possible)))
    assert not rule(_ResultView(_summary("FRAGILE", None, possible)))
    # flagged but the verdict was not downgraded
    assert not rule(_ResultView(_summary("ROBUST", "ROBUST", possible)))
    # not flagged
    none = {"mode": "B", "label": "NO TIME CONFOUND DETECTED"}
    assert not rule(_ResultView(_summary("FRAGILE", "ROBUST", none)))
    assert not rule(_ResultView(_summary("FRAGILE", "ROBUST", None)))


def test_mode_a_rule_needs_held_unadjusted_effect_removed_and_not_supported(kang_small):
    rule = _time_cases(plant_time_signal(kang_small, seed=0))["time_mode_a"].passes
    removed = {
        "mode": "A",
        "label": "POSSIBLE TIME CONFOUND",
        "effect_removed": True,
        "held_unadjusted": True,
    }
    assert rule(_ResultView(_summary("NOT SUPPORTED", "ROBUST", removed)))
    not_held = dict(removed, held_unadjusted=False)
    assert not rule(_ResultView(_summary("NOT SUPPORTED", "ROBUST", not_held)))
    assert not rule(
        _ResultView(_summary("NOT SUPPORTED", "ROBUST", dict(removed, held_unadjusted=None)))
    )
    assert not rule(
        _ResultView(_summary("NOT SUPPORTED", "ROBUST", dict(removed, effect_removed=False)))
    )
    assert not rule(_ResultView(_summary("FRAGILE", "ROBUST", removed)))
    assert not rule(_ResultView(_summary("NOT SUPPORTED", "ROBUST", None)))


def test_rule_texts_state_the_mechanical_conditions(kang_small):
    cases = _time_cases(plant_time_signal(kang_small, seed=0))
    b, a = cases["time_mode_b"], cases["time_mode_a"]
    assert "label_before_time_check is ROBUST" in b.passes_rule
    assert "POSSIBLE TIME CONFOUND" in b.passes_rule and "not ROBUST" in b.passes_rule
    assert "held_unadjusted True" in a.passes_rule and "effect_removed True" in a.passes_rule
    assert a.passes_rule.endswith("the verdict is NOT SUPPORTED")


def test_doublet_rule_requires_the_claim_to_hold_somewhere(sim):
    ref = sim.obs["cell_type"].astype(str)
    case = plant_doublets(sim, ref, type_a="Mono", type_b="B", n=5, seed=0)[0]

    def view(verdict: str, runs: list[tuple[bool, bool]]) -> _ResultView:
        return _ResultView(
            {
                "verdict": verdict,
                "top_breakers": [],
                "runs": [{"holds": h, "params": {"doublet_removal": dr}} for dr, h in runs],
            }
        )

    assert case.caught(view("FRAGILE", [(False, True), (True, False)]))
    assert not case.caught(view("NOT SUPPORTED", [(False, False), (True, False)]))  # never held
    assert not case.caught(view("FRAGILE", [(False, True), (True, True)]))
    assert not case.caught(view("ROBUST", [(False, True), (True, True)]))
    assert "held in at least one run" in case.caught_rule
    assert case.expected == "not ROBUST; holds only when doublet removal is off"


# ---------------------------------------------------------------- provenance


def test_time_signal_direction_matches_the_shipped_day_night_table():
    table = load_rhythmic_table()
    assert table is not None
    check = time_signal_day_night(table)
    assert set(check) == set(TIME_SIGNAL)  # every planted gene is on the S8 list
    disagree = sorted(g for g, r in check.items() if not r["agrees"])
    assert disagree == sorted(TIME_SIGNAL_S8_DISAGREE)
    for gene in ("PER1", "PER2", "PER3", "CRY2", "DBP", "TEF", "CIART"):
        assert check[gene]["expected"] == "day" and check[gene]["agrees"], gene
    for gene in ("ARNTL", "NPAS2", "CLOCK", "SQSTM1"):
        assert check[gene]["expected"] == "night" and check[gene]["agrees"], gene
    assert check["CRY1"] == {**check["CRY1"], "expected": "night", "observed": "day"}
    assert check["NR1D1"] == {**check["NR1D1"], "expected": "day", "observed": "night"}
    assert check["NR1D2"]["expected"] == "day" and check["NR1D2"]["observed"] != "day"
    assert time_signal_day_night(None) == {}


# -------------------------------------------------------------------- runner


def test_compact_time_check_keeps_what_the_rules_read():
    assert _compact_time_check(None) is None
    full = {
        "mode": "A",
        "label": "POSSIBLE TIME CONFOUND",
        "effect_removed": True,
        "effect_shrunk": False,
        "reasons": ["r"],
        "balance": {
            "p_permutation": 0.03,
            "circular_mean_difference_hours": -6.5,
            "unit_column": "replicate",
            "imbalanced": True,
            "imbalance_basis": "permutation",
            "times_by_sample": {"p1": 8.0},
        },
        "adjustment": {
            "held_unadjusted": True,
            "estimable": True,
            "survives_adjustment": False,
            "per_gene": [{"gene": "SQSTM1"}],
        },
        "panel": None,
        "rhythmic": {"available": True},
    }
    c = _compact_time_check(full)
    assert c is not None
    assert c["held_unadjusted"] is True and c["estimable"] is True and c["effect_removed"] is True
    assert c["imbalanced"] is True and c["imbalance_basis"] == "permutation"
    assert c["balance_unit_column"] == "replicate" and c["panel_unit_column"] is None
    assert c["balance_p"] == 0.03 and c["adjustment_per_gene"] == [{"gene": "SQSTM1"}]
    assert "times_by_sample" not in c and "rhythmic" not in c
    mode_b = _compact_time_check(
        {
            "mode": "B",
            "label": "NO TIME CONFOUND DETECTED",
            "panel": {
                "genes_usable": ["PER1", "ARNTL"],
                "n_panel_genes_different": 0,
                "unit_column": "replicate",
            },
        }
    )
    assert mode_b is not None
    assert mode_b["panel_genes_usable"] == 2 and mode_b["panel_unit_column"] == "replicate"
    assert mode_b["held_unadjusted"] is None and mode_b["estimable"] is None


def test_results_table_marks_planted_time_rows():
    def row(passes: bool, label: str) -> dict[str, Any]:
        return {
            "title": "time",
            "claim": "c",
            "expected": "e",
            "verdict": "FRAGILE",
            "n_holds": 12,
            "n_runs": 20,
            "cell_level_survival": 0.8,
            "pseudoreplication_flag": False,
            "passes": passes,
            "role": "planted_time",
            "part": "kang",
            "claim_type": "differential_expression",
            "time_check": {"label": label},
        }

    cases = [row(True, "POSSIBLE TIME CONFOUND"), row(False, "NOT ASSESSABLE")]
    lines = results_table({"cases": cases}).splitlines()[2:]
    assert lines[0].endswith("| 80% | POSSIBLE | caught |")
    assert lines[1].endswith("| 80% | n/a | MISSED |")


def test_kang_jobs_carry_their_part(kang_small, tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    kang_small.write_h5ad(data_dir / "kang_2018.h5ad")
    monkeypatch.setattr(runner, "download_kang", lambda path: path)
    monkeypatch.setattr(runner, "doublet_check", lambda geo_dir: {"skipped": "no GEO data"})
    jobs, meta, info = runner.build_kang_jobs(
        data_dir, tmp_path / "out", GridConfig(), Thresholds()
    )
    assert {j["part"] for j in jobs} == {"kang"}
    assert [j["key"] for j in jobs][-2:] == ["time_mode_b", "time_mode_a"]
    time_jobs = [j for j in jobs if j["key"].startswith("time_mode")]
    assert all(j["spec"]["columns"]["sample"] == ["replicate"] for j in time_jobs)
    assert meta["time_mode_a"]["role"] == "planted_time"
    assert set(info["planted_time"]) == {"groups", "collection_time_by_donor", "genes"}
    assert (data_dir / "kang_time_planted.h5ad").exists()


def test_every_benchmark_claim_has_its_own_report_folder(tmp_path):
    from recluster.engine import slugify

    cases = kang_cases(
        tmp_path / "real.h5ad",
        tmp_path / "split.h5ad",
        tmp_path / "shuffle.h5ad",
        {"A": ["d0"], "B": ["d1"]},
        tmp_path / "time.h5ad",
        {"d0": "A", "d1": "B"},
    )
    slugs = [slugify(c.spec.name) for c in cases]
    assert len(set(slugs)) == len(slugs), slugs
