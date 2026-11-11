"""The rules added after the first benchmark: QC cap, within-sample matching, per-marker
majority, label-or-marker candidates, pairing, the both-groups precondition, Scrublet fallback."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from recluster.config import GridConfig, PipelineParams, Thresholds
from recluster.data import detect_mito_prefix, extract_counts
from recluster.evaluators.composition import evaluate_composition, group_presence
from recluster.evaluators.pairing import pair_key
from recluster.matching import (
    lognorm_matrix,
    marker_score,
    match_population,
    per_gene_aurocs,
    within_sample_auroc,
)
from recluster.pipeline import Pipeline
from recluster.spec import Comparison, make_sample_key

from .conftest import mono_markers

# ------------------------------------------------------------------- QC cap


def test_qc_cap_pulls_back_a_threshold_that_removes_too_many_cells(sim):
    counts, _ = extract_counts(sim)
    pipe = Pipeline(counts, GridConfig(), mito_prefix=detect_mito_prefix(counts.var_names))
    n_genes = counts.obs["n_genes_by_counts"].to_numpy()
    harsh = int(np.quantile(n_genes, 0.6))  # would remove ~60 % of cells
    mg, mito, capped = pipe.effective_qc(PipelineParams(min_genes=harsh))
    assert "min_genes" in capped and mg < harsh
    assert abs((n_genes < mg).mean() - 0.25) < 0.03
    mg2, _, capped2 = pipe.effective_qc(PipelineParams(min_genes=50))
    assert mg2 == 50 and capped2 == {}
    stage = pipe.qc_stage(PipelineParams(min_genes=harsh))
    assert stage.effective["capped"] == capped
    assert any("capped" in w for w in stage.warnings)
    run = pipe.run(PipelineParams(min_genes=harsh, batch_correction=False))
    assert run.effective["min_genes"] == mg


def test_qc_cap_is_configurable(sim):
    counts, _ = extract_counts(sim)
    pipe = Pipeline(counts, GridConfig(qc_max_cell_loss=1.0), mito_prefix="MT-")
    n_genes = counts.obs["n_genes_by_counts"].to_numpy()
    harsh = int(np.quantile(n_genes, 0.6))
    assert pipe.effective_qc(PipelineParams(min_genes=harsh))[2] == {}


# ----------------------------------------------------------- matching rules


def test_within_sample_auroc_ignores_a_sample_level_shift():
    rng = np.random.default_rng(0)
    # cluster "pop" scores high in every sample; one whole sample is shifted up, as a QC
    # change that keeps more low-quality cells in that sample can do
    samples = np.repeat(["s1", "s2"], 200)
    groups = np.array(["pop"] * 50 + ["rest"] * 150 + ["pop"] * 50 + ["rest"] * 150)
    score = np.where(groups == "pop", rng.normal(3, 1, 400), rng.normal(0, 1, 400))
    before = within_sample_auroc(score, groups, samples)["pop"]
    score[200:] += 2.5
    after = within_sample_auroc(score, groups, samples)["pop"]
    assert before > 0.9 and abs(after - before) < 1e-9
    one_sample = within_sample_auroc(score[:200], groups[:200], samples[:200])
    assert "pop" in one_sample and "rest" in one_sample


def test_per_marker_majority_rejects_a_set_carried_by_half_its_genes(sim, truth):
    counts, _ = extract_counts(sim)
    clusters = sim.obs["cell_type"].astype(str).rename("cluster")
    mono, b = mono_markers(truth)[:2], list(truth["markers"]["B"][:2])
    frame = lognorm_matrix(counts, mono + b)
    score = marker_score(frame, mono + b, [])
    thr = Thresholds()
    m = match_population(
        clusters,
        score=score,
        reference=None,
        sample=None,
        min_cells=20,
        thr=thr,
        marker_frame=frame,
    )
    assert m.clusters == []  # neither Mono nor B carries more than half of the four markers
    pg = per_gene_aurocs(frame.loc[clusters.index], clusters.to_numpy())
    assert (pg.loc["Mono", mono] > 0.8).all() and (pg.loc["Mono", b] < 0.6).all()
    m_all = match_population(
        clusters,
        score=marker_score(frame, mono + b, [], logic="all"),
        reference=None,
        sample=None,
        min_cells=20,
        thr=thr,
        marker_frame=frame,
        logic="all",
    )
    assert m_all.clusters == []


def test_candidate_rule_either_recovers_a_labelled_population_with_weak_markers(sim, truth):
    counts, _ = extract_counts(sim)
    clusters = sim.obs["cell_type"].astype(str).rename("cluster")
    ref = clusters.index[clusters == "Mono"]
    weak = [truth["null_genes"][0]]  # a gene with no enrichment at all
    frame = lognorm_matrix(counts, weak)
    score = marker_score(frame, weak, [])
    thr = Thresholds()
    either = match_population(
        clusters,
        score=score,
        reference=ref,
        sample=None,
        min_cells=20,
        thr=thr,
        marker_frame=frame,
        reference_is_labelled=True,
    )
    assert either.clusters == ["Mono"] and either.candidate_by == {"Mono": "labels"}
    markers_only = match_population(
        clusters,
        score=score,
        reference=ref,
        sample=None,
        min_cells=20,
        thr=Thresholds(candidate_rule="markers"),
        marker_frame=frame,
        reference_is_labelled=True,
    )
    assert markers_only.clusters == []
    # a reference derived from an earlier run never makes a cluster a candidate by itself
    derived = match_population(
        clusters, score=score, reference=ref, sample=None, min_cells=20, thr=thr, marker_frame=frame
    )
    assert derived.clusters == [] and derived.mode == "markers"


def test_jaccard_on_kept_cells_reports_qc_removal(sim, truth):
    counts, _ = extract_counts(sim)
    clusters = sim.obs["cell_type"].astype(str).rename("cluster")
    ref = clusters.index[clusters == "Mono"]
    kept = clusters.iloc[: len(clusters) // 2]  # pretend QC removed half of the cells
    score = marker_score(lognorm_matrix(counts, mono_markers(truth)), mono_markers(truth), [])
    m = match_population(
        kept, score=score, reference=ref, sample=None, min_cells=20, thr=Thresholds()
    )
    assert m.jaccard is not None and m.jaccard > 0.95
    assert m.reference_qc_removed is not None and 0.3 < m.reference_qc_removed < 0.7
    old = match_population(
        kept,
        score=score,
        reference=ref,
        sample=None,
        min_cells=20,
        thr=Thresholds(jaccard_on_kept_cells=False),
    )
    assert old.jaccard is not None and old.jaccard < 0.7


# -------------------------------------------------------------- pairing


def test_pairing_detected_only_when_every_donor_is_in_both_groups(sim):
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cond = sim.obs["condition"].astype(str)
    pair = pair_key(
        sample, cond, sim.obs_names, "stim", "ctrl", ["donor", "condition"], "condition", sim.obs
    )
    assert pair is not None and pair.nunique() == 6
    fake = (
        sim.obs["donor"]
        .astype(str)
        .isin(["donor0", "donor1", "donor2"])
        .map({True: "A", False: "B"})
    )
    obs = sim.obs.assign(fake=fake.to_numpy())
    none = pair_key(sample, obs["fake"], sim.obs_names, "A", "B", ["donor", "fake"], "fake", obs)
    assert none is None
    assert (
        pair_key(sample, cond, sim.obs_names, "stim", "ctrl", ["donor"], "condition", sim.obs)
        is None
    )


def test_paired_composition_uses_per_donor_differences(sim):
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cond = sim.obs["condition"].astype(str)
    pair = pair_key(
        sample, cond, sim.obs_names, "stim", "ctrl", ["donor", "condition"], "condition", sim.obs
    )
    mono = sim.obs_names[(sim.obs["cell_type"] == "Mono").to_numpy()]
    cmp = Comparison(group_a="stim", group_b="ctrl", direction="higher_in_a")
    paired = evaluate_composition(mono, sim.obs_names, sample, cond, cmp, Thresholds(), pair=pair)
    unpaired = evaluate_composition(mono, sim.obs_names, sample, cond, cmp, Thresholds())
    assert paired["design"] == "paired" and unpaired["design"] == "unpaired"
    assert paired["holds"] and unpaired["holds"]
    assert paired["p_logit_ttest"] != unpaired["p_logit_ttest"]


# -------------------------------------------- population in both groups


def test_composition_not_testable_when_population_missing_from_one_group(sim):
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cond = sim.obs["condition"].astype(str)
    mono_ctrl_only = sim.obs_names[
        ((sim.obs["cell_type"] == "Mono") & (sim.obs["condition"] == "ctrl")).to_numpy()
    ]
    cmp = Comparison(group_a="stim", group_b="ctrl", direction="different")
    res = evaluate_composition(mono_ctrl_only, sim.obs_names, sample, cond, cmp, Thresholds())
    assert res["evaluable"] is False and res["holds"] is False
    assert any("cannot test" in r for r in res["reasons"])
    table = pd.DataFrame(
        {
            "n": [100, 100, 100, 100],
            "k": [30, 25, 0, 2],
            "condition": ["stim", "stim", "ctrl", "ctrl"],
        },
        index=["a1", "a2", "b1", "b2"],
    )
    ok, reasons = group_presence(table, cmp, Thresholds())
    assert not ok and "ctrl" in reasons[0]


# ------------------------------------------------------- Scrublet fallback


def test_scrublet_fallback_calls_expected_fraction(monkeypatch, sim):
    import scanpy as sc

    counts, _ = extract_counts(sim)
    pipe = Pipeline(counts, GridConfig(), mito_prefix="MT-")

    def fake_scrublet(adata, **kw):  # a threshold above every score: Scrublet calls nothing
        rng = np.random.default_rng(0)
        adata.obs["doublet_score"] = rng.uniform(0, 0.3, adata.n_obs)
        adata.obs["predicted_doublet"] = False

    monkeypatch.setattr(sc.pp, "scrublet", fake_scrublet)
    idx = np.arange(counts.n_obs)
    scores, pred, info = pipe._scrublet_subset(counts, idx)
    assert info["fallback"] is True
    expected = min(0.30, GridConfig().scrublet_expected_rate_per_1000 * counts.n_obs / 1000)
    assert pred.sum() == round(expected * counts.n_obs)
    assert scores[pred].min() >= scores[~pred].max()  # the top-scoring cells were called


def test_scrublet_fallback_not_used_when_calls_are_normal(monkeypatch, sim):
    import scanpy as sc

    counts, _ = extract_counts(sim)
    pipe = Pipeline(counts, GridConfig(), mito_prefix="MT-")

    def fake_scrublet(adata, **kw):
        rng = np.random.default_rng(1)
        adata.obs["doublet_score"] = rng.uniform(0, 0.3, adata.n_obs)
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > 0.29  # ~3 % called

    monkeypatch.setattr(sc.pp, "scrublet", fake_scrublet)
    _, pred, info = pipe._scrublet_subset(counts, np.arange(counts.n_obs))
    assert info["fallback"] is False and pred.sum() > 0


@pytest.mark.parametrize("logic", ["any", "all"])
def test_marker_logic_accepted_by_spec(logic):
    from recluster.spec import ClaimSpec

    spec = ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": "x",
            "target": {"markers": ["A", "B"], "marker_logic": logic},
        }
    )
    assert spec.target.marker_logic == logic
