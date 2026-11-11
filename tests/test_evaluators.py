"""Each evaluator: a claim that is true by construction passes, a false one fails.

These tests bypass clustering and hand the evaluators the planted truth labels, so they test
the statistics, not the pipeline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from recluster.data import extract_counts
from recluster.evaluators.composition import evaluate_composition
from recluster.evaluators.de import bh_adjust, evaluate_de
from recluster.evaluators.population import evaluate_population
from recluster.matching import lognorm_matrix, marker_score, match_population
from recluster.spec import Comparison, make_sample_key

from .conftest import mono_markers


@pytest.fixture(scope="module")
def ctx(sim):
    counts, _ = extract_counts(sim)
    truth = dict(sim.uns["truth"])
    genes = mono_markers(truth) + list(truth["de_genes"]) + list(truth["null_genes"])
    genes = list(dict.fromkeys(genes))
    return {
        "counts": counts,
        "truth": truth,
        "sample": make_sample_key(sim.obs, ["donor", "condition"]),
        "condition": sim.obs["condition"].astype(str),
        "labels": sim.obs["cell_type"].astype(str),
        "lognorm": lognorm_matrix(counts, genes),
    }


# ------------------------------------------------------------- population


def test_population_recovered_from_truth_clusters(ctx, thresholds):
    clusters = ctx["labels"].rename("cluster")  # perfect clustering
    score = marker_score(ctx["lognorm"], mono_markers(ctx["truth"]), [])
    ref = ctx["labels"].index[ctx["labels"] == "Mono"]
    m = match_population(
        clusters, score=score, reference=ref, sample=ctx["sample"], min_cells=20, thr=thresholds
    )
    assert m.clusters == ["Mono"]
    assert m.jaccard == 1.0 and m.marker_auroc > 0.95
    assert m.n_samples_present == 12
    ev = evaluate_population(m, thresholds)
    assert ev["holds"] and ev["recovered"] and ev["replicate_ok"]


def test_population_fails_for_random_markers(ctx, thresholds):
    clusters = ctx["labels"].rename("cluster")
    rng = np.random.default_rng(0)
    pick = rng.choice(ctx["counts"].var_names, 8, replace=False)
    fake = [g for g in pick if g.startswith("G")][:5]
    score = marker_score(lognorm_matrix(ctx["counts"], fake), fake, [])
    m = match_population(
        clusters, score=score, reference=None, sample=ctx["sample"], min_cells=20, thr=thresholds
    )
    assert not m.found
    assert m.best_cluster_auroc < thresholds.marker_auroc_min
    assert evaluate_population(m, thresholds)["holds"] is False


def test_population_fails_when_confined_to_one_sample(ctx, thresholds):
    clusters = ctx["labels"].rename("cluster")
    score = marker_score(ctx["lognorm"], mono_markers(ctx["truth"]), [])
    one_sample = ctx["sample"].iloc[0]
    sample = ctx["sample"].copy()
    mono = ctx["labels"] == "Mono"
    sample[mono.to_numpy()] = one_sample  # pretend every Mono cell came from one sample
    m = match_population(
        clusters, score=score, reference=None, sample=sample, min_cells=20, thr=thresholds
    )
    assert m.found and m.n_samples_present == 1
    ev = evaluate_population(m, thresholds)
    assert ev["recovered"] and not ev["replicate_ok"] and not ev["holds"]


def test_label_mode_matching_tolerates_split_clusters(ctx, thresholds):
    labels = ctx["labels"].copy()
    mono_idx = labels.index[labels == "Mono"]
    half = mono_idx[: len(mono_idx) // 2]
    clusters = labels.copy()
    clusters.loc[half] = "Mono_a"
    clusters.loc[mono_idx.difference(half)] = "Mono_b"
    m = match_population(
        clusters.rename("cluster"),
        score=None,
        reference=mono_idx,
        sample=ctx["sample"],
        min_cells=20,
        thr=thresholds,
    )
    assert sorted(m.clusters) == ["Mono_a", "Mono_b"]
    assert m.jaccard == 1.0


# ------------------------------------------------------------ composition


def test_composition_detects_planted_shift(ctx, thresholds):
    mono = ctx["labels"].index[ctx["labels"] == "Mono"]
    res = evaluate_composition(
        mono,
        ctx["labels"].index,
        ctx["sample"],
        ctx["condition"],
        Comparison(group_a="stim", group_b="ctrl", direction="higher_in_a"),
        thresholds,
    )
    assert res["holds"] and res["direction_ok"]
    assert res["p_logit_ttest"] < 0.01
    assert res["effect_log_odds_a_minus_b"] > 0
    assert res["n_samples_a"] == 6 and res["n_samples_b"] == 6


def test_composition_rejects_wrong_direction(ctx, thresholds):
    mono = ctx["labels"].index[ctx["labels"] == "Mono"]
    res = evaluate_composition(
        mono,
        ctx["labels"].index,
        ctx["sample"],
        ctx["condition"],
        Comparison(group_a="stim", group_b="ctrl", direction="lower_in_a"),
        thresholds,
    )
    assert not res["holds"]


def test_composition_no_shift_fails_replicate_aware(ctx, thresholds):
    t_cells = ctx["labels"].index[ctx["labels"] == "T"]
    res = evaluate_composition(
        t_cells,
        ctx["labels"].index,
        ctx["sample"],
        ctx["condition"],
        Comparison(group_a="stim", group_b="ctrl", direction="higher_in_a"),
        thresholds,
    )
    # T shrinks only as a side effect of Mono growing; direction asserted here is wrong anyway
    assert not res["holds"]


def test_composition_too_few_replicates_is_undefined(ctx, thresholds):
    mono = ctx["labels"].index[ctx["labels"] == "Mono"]
    one_donor = ctx["sample"].str.startswith("donor0")
    cells = ctx["labels"].index[one_donor.to_numpy()]
    res = evaluate_composition(
        mono.intersection(cells),
        cells,
        ctx["sample"],
        ctx["condition"],
        Comparison(group_a="stim", group_b="ctrl", direction="higher_in_a"),
        thresholds,
    )
    assert res["p_logit_ttest"] is None and not res["holds"]
    assert any("fewer than 2 replicates" in w for w in res["warnings"])


# --------------------------------------------------------------------- DE


def _de(ctx, thresholds, genes, direction="up_in_a", require="all"):
    mono = ctx["labels"].index[ctx["labels"] == "Mono"]
    return evaluate_de(
        cells=mono,
        counts=ctx["counts"],
        lognorm=ctx["lognorm"],
        sample=ctx["sample"],
        condition=ctx["condition"],
        comparison=Comparison(group_a="stim", group_b="ctrl", direction=direction),
        genes=genes,
        require=require,
        thr=thresholds,
        n_cpus=2,
    )


def test_de_planted_genes_pass_both_levels(ctx, thresholds):
    res = _de(ctx, thresholds, list(ctx["truth"]["de_genes"]))
    assert res["holds"] and res["holds_cell_level"]
    per = pd.DataFrame(res["per_gene"])
    assert (per["agreement"] == "both").all()
    assert (per["pb_log2fc"] > 1).all()
    assert res["n_samples_a"] == 6 and res["n_samples_b"] == 6


def test_de_null_genes_fail(ctx, thresholds):
    res = _de(ctx, thresholds, list(ctx["truth"]["null_genes"]), direction="either")
    assert not res["holds"]
    per = pd.DataFrame(res["per_gene"])
    assert (per["pb_padj"] > thresholds.alpha).all()


def test_de_direction_matters(ctx, thresholds):
    res = _de(ctx, thresholds, list(ctx["truth"]["de_genes"]), direction="down_in_a")
    assert not res["holds"] and not res["holds_cell_level"]


def test_de_require_any_vs_all(ctx, thresholds):
    genes = [ctx["truth"]["de_genes"][0], ctx["truth"]["null_genes"][0]]
    assert not _de(ctx, thresholds, genes, direction="either", require="all")["holds"]
    assert _de(ctx, thresholds, genes, direction="either", require="any")["holds"]


def test_bh_adjust_matches_known_values():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan])
    adj = bh_adjust(p)
    assert np.isnan(adj[4])
    # sorted p: .01,.03,.04,.2 -> .04,.06,.0533,.2 -> running min from the right
    expected = np.array([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.2])
    assert np.allclose(adj[:4], expected)
