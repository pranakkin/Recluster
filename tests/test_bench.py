"""Benchmark constructions: deterministic, correctly shaped, and honest in their bookkeeping.

The real-data runs are not exercised here (they need downloads and minutes of compute);
these tests cover the pieces that must be right for the benchmark to mean anything.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from recluster.bench.doublets import _agreement, _unique
from recluster.bench.kang import check_kang, donor_level_shuffle, fake_split_within_control
from recluster.bench.runner import _ResultView, results_table


@pytest.fixture
def kang_like(sim) -> AnnData:
    """The synthetic dataset with Kang's column names."""
    a = sim.copy()
    a.obs = a.obs.rename(columns={"condition": "label", "donor": "replicate"})
    a.obs["replicate"] = ("patient_" + a.obs["replicate"].astype(str)).astype("category")
    return a


def test_check_kang_reports_columns(kang_like):
    info = check_kang(kang_like)
    assert info["conditions"] == ["ctrl", "stim"] and len(info["donors"]) == 6
    assert info["doublet_columns"] == []
    assert info["n_mito_genes"] == 10


def test_check_kang_refuses_wrong_columns(kang_like):
    bad = kang_like.copy()
    bad.obs = bad.obs.drop(columns=["replicate"])
    with pytest.raises(RuntimeError, match="replicate"):
        check_kang(bad)


def test_fake_split_is_control_only_and_balanced(kang_like):
    split = fake_split_within_control(kang_like, seed=0)
    assert set(split.obs["label"].astype(str)) == {"ctrl"}
    groups = split.uns["fake_split"]
    assert len(groups["A"]) == 3 and len(groups["B"]) == 3
    assert not set(groups["A"]) & set(groups["B"])
    per_donor = split.obs.groupby("replicate", observed=True)["fake_group"].nunique()
    assert (per_donor == 1).all()  # a donor is in exactly one fake group
    again = fake_split_within_control(kang_like, seed=0)
    assert again.uns["fake_split"] == groups
    assert fake_split_within_control(kang_like, seed=1).uns["fake_split"] != groups


def test_donor_shuffle_is_balanced_by_condition(kang_like):
    sh = donor_level_shuffle(kang_like, seed=0)
    tab = sh.obs.groupby(["fake_group", "label"], observed=True)["replicate"].nunique().unstack()
    assert (tab == 3).all().all()  # each fake group: 3 ctrl + 3 stim samples
    key = sh.obs["replicate"].astype(str) + "|" + sh.obs["label"].astype(str)
    assert (sh.obs.groupby(key.to_numpy())["fake_group"].nunique() == 1).all()


def test_doublet_agreement_metrics():
    obs = pd.DataFrame(
        {
            "demuxlet": ["doublet", "doublet", "singlet", "singlet", "ambs", "singlet"],
            "predicted_doublet": [True, False, True, False, True, False],
            "doublet_score": [0.9, 0.4, 0.5, 0.1, 0.8, 0.2],
        }
    )
    a = _agreement(obs)
    assert a["n_cells"] == 5 and a["ambiguous_excluded"] == 1
    assert (a["tp"], a["fp"], a["fn"], a["tn"]) == (1, 1, 1, 2)
    assert a["precision"] == 0.5 and a["recall"] == 0.5
    assert 0 < a["auroc_of_score"] <= 1


def test_unique_gene_names():
    assert _unique(pd.Index(["A", "B", "A", "A"])).tolist() == ["A", "B", "A-1", "A-2"]


def test_results_table_marks_misses_and_false_alarms():
    results = {
        "cases": [
            {
                "title": "t",
                "claim": "c",
                "expected": "ROBUST",
                "verdict": "FRAGILE",
                "n_holds": 10,
                "n_runs": 20,
                "cell_level_survival": 0.9,
                "pseudoreplication_flag": True,
                "passes": False,
                "role": "true",
                "part": "kang",
            },
            {
                "title": "p",
                "claim": "c",
                "expected": "not ROBUST",
                "verdict": "ROBUST",
                "n_holds": 20,
                "n_runs": 20,
                "cell_level_survival": None,
                "pseudoreplication_flag": False,
                "passes": False,
                "role": "planted",
                "part": "planted",
            },
        ]
    }
    md = results_table(results)
    assert "MISSED" in md and "⚑" in md and "90%" in md


def test_result_view_exposes_rule_inputs():
    s = {
        "verdict": "FRAGILE",
        "top_breakers": [{"axis": "resolution"}],
        "runs": [
            {"holds": True, "params": {"doublet_removal": False}},
            {"holds": False, "params": {"doublet_removal": True}},
        ],
    }
    v = _ResultView(s)
    assert v.verdict.label == "FRAGILE" and v.verdict.breakers[0]["axis"] == "resolution"
    assert [r.holds for r in v.runs] == [True, False]


PBMC3K = Path("data/pbmc3k_raw.h5ad")


@pytest.mark.skipif(not PBMC3K.exists(), reason="PBMC3k not downloaded")
def test_planted_constructions_are_deterministic_and_labelled():
    from recluster.bench.pbmc3k import load_pbmc3k, reference_labels
    from recluster.bench.planted import plant_doublets, plant_qc_deletion

    a = load_pbmc3k("data")
    ref = reference_labels(a)
    assert {"CD4 T", "CD8 T", "NK", "B", "CD14 Mono"} <= set(ref.unique())
    d1, d2 = plant_doublets(a, ref, n=20, seed=0)[0], plant_doublets(a, ref, n=20, seed=0)[0]
    assert d1.adata.n_obs == a.n_obs + 20
    assert np.array_equal(d1.adata.X[-20:].toarray(), d2.adata.X[-20:].toarray())
    assert d1.adata.obs["injected_doublet"].sum() == 20
    q = plant_qc_deletion(a, ref, seed=0)[0]
    assert q.adata.n_obs == a.n_obs and q.caught_rule.startswith("top entry")
