"""Matching rules: the stratified within-sample AUROC, which statistic judged each candidate, and
when a reference set may make a cluster a candidate by labels."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm, rankdata

from recluster.config import Thresholds
from recluster.data import extract_counts
from recluster.matching import (
    Match,
    auroc_by_group,
    lognorm_matrix,
    marker_score,
    match_population,
    within_sample_auroc,
)

from .conftest import mono_markers

# a rare population: N_POP cells in each of N_SAMPLES samples, next to N_REST other cells
N_SAMPLES, N_POP, N_REST = 16, 12, 300
TRUE_AUROC = 0.85
MU = float(np.sqrt(2) * norm.ppf(TRUE_AUROC))  # normal(MU, 1) vs normal(0, 1) has this AUROC


def _rare_population(rng):
    per = N_POP + N_REST
    samples = np.repeat([f"s{i:02d}" for i in range(N_SAMPLES)], per)
    groups = np.tile(np.array(["pop"] * N_POP + ["rest"] * N_REST), N_SAMPLES)
    score = rng.normal(0.0, 1.0, len(groups)) + np.where(groups == "pop", MU, 0.0)
    return score, groups, samples


def _min_over_samples(score, groups, samples, cluster="pop"):
    """The statistic the stratified AUROC replaced: the worst per-sample AUROC of the cluster."""
    worst = 1.0
    for s in np.unique(samples):
        in_s = samples == s
        ranks = rankdata(score[in_s])
        mask = groups[in_s] == cluster
        n1, n0 = int(mask.sum()), int((~mask).sum())
        worst = min(worst, (ranks[mask].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))
    return worst


def _as_series(score, groups, samples):
    index = pd.Index([f"cell{i}" for i in range(len(groups))])
    return (
        pd.Series(groups, index=index, name="cluster"),
        pd.Series(score, index=index),
        pd.Series(samples, index=index),
    )


# ------------------------------------------------ stratified within-sample AUROC


def test_stratified_auroc_pools_cell_pairs_not_samples():
    score = np.array([3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    groups = np.array(["pop", "pop", "rest", "pop", "rest", "rest", "rest", "rest", "only", "only"])
    samples = np.array(["a", "a", "a", "b", "b", "b", "b", "b", "c", "c"])
    out = within_sample_auroc(score, groups, samples)
    # sample a: U = 2 of 2 pairs; sample b: U = 0 of 4 pairs -> 2/6, not the mean of 1.0 and 0.0
    assert out["pop"] == pytest.approx(2 / 6)
    assert out["rest"] == pytest.approx(4 / 6)
    assert "only" not in out  # fills its sample, so nothing there to rank it against


def test_stratified_auroc_is_unbiased_and_keeps_a_rare_population():
    cutoff = Thresholds().marker_auroc_min
    stratified, worst = [], []
    for seed in range(200):
        score, groups, samples = _rare_population(np.random.default_rng(seed))
        stratified.append(within_sample_auroc(score, groups, samples)["pop"])
        worst.append(_min_over_samples(score, groups, samples))
    strat, mins = np.array(stratified), np.array(worst)
    assert abs(strat.mean() - TRUE_AUROC) < 0.03
    assert (strat < cutoff).mean() < 0.10
    assert (mins < cutoff).mean() > 0.50  # the worst of 16 samples of 12 cells is usually low


def test_stratified_auroc_ignores_a_sample_level_shift():
    score, groups, samples = _rare_population(np.random.default_rng(0))
    base = within_sample_auroc(score, groups, samples)["pop"]
    shifted = score.copy()
    # every score of one sample sits higher, as a batch effect does
    shifted[samples == "s03"] += 2.5
    after = within_sample_auroc(shifted, groups, samples)["pop"]
    assert after == pytest.approx(base)  # ranks inside the sample did not change
    # the whole-dataset AUROC moves: that sample's other cells now outrank the other populations
    assert auroc_by_group(shifted, groups)["pop"] < auroc_by_group(score, groups)["pop"] - 0.02


# --------------------------------------------------- what judged a candidate


def test_judged_by_and_marker_basis_name_the_statistic():
    score, groups, samples = _rare_population(np.random.default_rng(1))
    clusters, s, sample = _as_series(score, groups, samples)
    thr = Thresholds()
    within = match_population(
        clusters, score=s, reference=None, sample=sample, min_cells=20, thr=thr
    )
    assert within.clusters == ["pop"] and within.mode == "markers"
    assert within.judged_by == {"pop": "within-sample"}
    assert within.marker_basis == "within-sample"
    expected = within_sample_auroc(score, groups, samples)["pop"]
    assert within.best_cluster_auroc == pytest.approx(expected)
    overall = match_population(
        clusters, score=s, reference=None, sample=None, min_cells=20, thr=thr
    )
    assert overall.judged_by == {"pop": "overall"} and overall.marker_basis == "overall"
    assert overall.best_cluster_auroc == pytest.approx(auroc_by_group(score, groups)["pop"])


def test_marker_basis_is_mixed_when_a_candidate_fills_its_samples():
    score, groups, samples = _rare_population(np.random.default_rng(2))
    # a second population that is the only content of its own two samples: no sample can judge
    # it, so the whole-dataset AUROC does
    extra = 40
    score = np.concatenate([score, np.random.default_rng(3).normal(4.0, 1.0, extra)])
    groups = np.concatenate([groups, np.repeat("solo", extra)])
    samples = np.concatenate([samples, np.repeat(["t0", "t1"], extra // 2)])
    clusters, s, sample = _as_series(score, groups, samples)
    m = match_population(
        clusters, score=s, reference=None, sample=sample, min_cells=20, thr=Thresholds()
    )
    assert m.clusters == ["pop", "solo"]
    assert m.judged_by == {"pop": "within-sample", "solo": "overall"}
    assert m.marker_basis == "mixed"


def test_marker_basis_without_candidates_names_the_statistic_behind_the_best_value():
    rng = np.random.default_rng(4)
    _, groups, samples = _rare_population(rng)
    score = rng.normal(0.0, 1.0, len(groups))  # no enrichment anywhere
    clusters, s, sample = _as_series(score, groups, samples)
    m = match_population(
        clusters, score=s, reference=None, sample=sample, min_cells=20, thr=Thresholds()
    )
    assert m.clusters == [] and m.judged_by == {} and m.marker_basis == "within-sample"
    assert m.best_cluster_auroc is not None and m.best_cluster_auroc < 0.8


def test_labels_only_match_has_no_marker_basis():
    score, groups, samples = _rare_population(np.random.default_rng(5))
    clusters, _, sample = _as_series(score, groups, samples)
    ref = clusters.index[clusters == "pop"]
    m = match_population(
        clusters,
        score=None,
        reference=ref,
        sample=sample,
        min_cells=20,
        thr=Thresholds(),
        reference_is_labelled=True,
    )
    assert m.mode == "labels" and m.clusters == ["pop"] and m.candidate_by == {"pop": "labels"}
    assert m.marker_basis is None and m.judged_by == {} and m.best_cluster_auroc is None
    # no markers and a reference nobody labelled: nothing to match by
    with pytest.raises(ValueError):
        match_population(
            clusters, score=None, reference=ref, sample=sample, min_cells=20, thr=Thresholds()
        )


# ------------------------------------------------------ reference provenance


def test_only_a_labelled_reference_makes_a_merged_cluster_a_candidate(sim, truth):
    counts, _ = extract_counts(sim)
    types = sim.obs["cell_type"].astype(str)
    sample = sim.obs["sample"].astype(str)
    rng = np.random.default_rng(0)
    mono = types.index[types == "Mono"]
    merged = pd.Index(rng.choice(mono, 40, replace=False)).append(
        pd.Index(rng.choice(types.index[types == "B"], 30, replace=False))
    )
    clusters = types.rename("cluster").copy()
    # a low-resolution cluster that swallowed part of the population
    clusters.loc[merged] = "merged"
    markers = mono_markers(truth)
    frame = lognorm_matrix(counts, markers)
    score = marker_score(frame, markers, [])
    thr = Thresholds()
    judged = within_sample_auroc(
        score.loc[clusters.index].to_numpy(), clusters.to_numpy(), sample.to_numpy()
    )
    assert judged["merged"] < thr.marker_auroc_min  # the markers alone do not support it
    carried = match_population(  # the reference came from an earlier run's match
        clusters,
        score=score,
        reference=mono,
        sample=sample,
        min_cells=20,
        thr=thr,
        marker_frame=frame,
    )
    assert carried.mode == "markers" and carried.clusters == ["Mono"]
    assert carried.candidate_by == {"Mono": "markers"}
    assert carried.jaccard is not None and carried.n_reference == len(mono)  # still compared with
    labelled = match_population(  # the reference came from the claim's label column
        clusters,
        score=score,
        reference=mono,
        sample=sample,
        min_cells=20,
        thr=thr,
        marker_frame=frame,
        reference_is_labelled=True,
    )
    assert labelled.mode == "either" and labelled.clusters == ["Mono", "merged"]
    assert labelled.candidate_by == {"Mono": "both", "merged": "labels"}
    assert labelled.judged_by == {"Mono": "within-sample"}


# ---------------------------------------------------------------- Match


def test_match_constructs_without_cells():
    m = Match(
        mode="markers",
        clusters=[],
        n_cells=0,
        found=False,
        best_cluster_auroc=None,
        marker_basis=None,
        marker_auroc=None,
        per_marker_auroc=None,
        markers_enriched=None,
        jaccard=None,
        precision=None,
        recall=None,
        n_reference=None,
        reference_qc_removed=None,
        n_samples_present=None,
        max_sample_share=None,
    )
    assert len(m.cells) == 0 and m.candidate_by == {} and m.judged_by == {}
    assert "cells" not in m.to_dict() and m.to_dict()["judged_by"] == {}
