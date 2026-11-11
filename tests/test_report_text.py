"""Report text: design-dependent labels, methods wording, and the informational blocks."""

from __future__ import annotations

import pytest

from recluster.config import GridConfig, Thresholds
from recluster.engine import run_claim
from recluster.report.html import _composition_table, _time_html, render_report
from recluster.report.methods import methods_paragraph
from recluster.results import ClaimResult, Verdict
from recluster.spec import ClaimSpec

from .conftest import mono_markers

COLS = {
    "condition": "condition",
    "sample": ["donor", "condition"],
    "batch": "donor",
    "labels": "cell_type",
}
COMPOSITION = {"group_a": "stim", "group_b": "ctrl", "direction": "higher_in_a"}
DE = {"group_a": "stim", "group_b": "ctrl", "direction": "up_in_a"}


def _spec(claim_type, truth, **extra):
    d = {
        "claim_type": claim_type,
        "name": f"{claim_type} text",
        "columns": COLS,
        "target": {"markers": mono_markers(truth), "label_value": "Mono"},
    }
    d.update(extra)
    return ClaimSpec.model_validate(d)


@pytest.fixture(scope="module")
def two_run_grid():
    """Baseline plus the first one-at-a-time variation."""
    return GridConfig(max_runs=2, n_random=0)


@pytest.fixture(scope="module")
def composition_result(sim_path, truth, two_run_grid):
    spec = _spec("composition_shift", truth, comparison=COMPOSITION)
    return run_claim(sim_path, spec, grid=two_run_grid, progress=None, write=False)


@pytest.fixture(scope="module")
def de_result(sim_path, truth, two_run_grid):
    spec = _spec("differential_expression", truth, comparison=DE, genes=list(truth["de_genes"][:2]))
    return run_claim(sim_path, spec, grid=two_run_grid, progress=None, write=False)


@pytest.fixture(scope="module")
def population_result(sim_path, truth, two_run_grid):
    spec = _spec("population_exists", truth)
    return run_claim(sim_path, spec, grid=two_run_grid, progress=None, write=False)


# ------------------------------------------------------------- composition labels


def test_composition_labels_follow_the_design():
    base = {
        "n_samples_a": 6,
        "n_samples_b": 6,
        "n_cells_a": 100,
        "n_cells_b": 100,
        "p_logit_ttest": 0.01,
        "p_mannwhitney": 0.02,
    }
    paired = _composition_table(dict(base, design="paired"), "x")
    assert "paired t on per-donor logit differences, p" in paired
    assert "Wilcoxon signed-rank on per-donor differences, p" in paired
    assert "Welch" not in paired and "Mann-Whitney" not in paired
    unpaired = _composition_table(dict(base, design="unpaired"), "x")
    assert "Welch t on logit proportions, p" in unpaired
    assert "Mann-Whitney U on proportions, p" in unpaired
    assert "paired t" not in unpaired and "signed-rank" not in unpaired


def test_composition_report_uses_paired_labels(composition_result):
    res = composition_result
    assert res.data_summary["paired_design"] is True
    assert res.runs[0].details["design"] == "paired"
    assert res.as_labelled["design"] == "paired"
    html = render_report(res)
    assert "paired t on per-donor logit differences, p" in html
    assert "Wilcoxon signed-rank on per-donor differences, p" in html
    assert "Welch t on logit proportions" not in html
    assert "Mann-Whitney U on proportions" not in html


# ------------------------------------------------------------------- methods text


def test_methods_scrublet_rate_is_a_percentage(composition_result):
    m = composition_result.methods
    assert "0.8% per 1,000 cells" in m
    assert "0.008 per" not in m


def test_methods_permutation_test_and_units(composition_result, de_result):
    permutation = (
        "a permutation test (exact when at most 20,000 label assignments exist, otherwise "
        "5,000 seeded draws)"
    )
    units = (
        "at the level of the sample, or of the donor when samples are nested in donors that "
        "share one condition and one collection time"
    )
    for res in (composition_result, de_result):
        assert permutation in res.methods
        assert units in res.methods
        assert "exact permutation test" not in res.methods


def test_methods_marker_match_is_stratified(population_result):
    m = population_result.methods
    assert "a stratified within-sample AUROC" in m
    assert "the Mann-Whitney U statistics of all samples pooled" in m
    assert "the minimum over samples" not in m


def test_methods_describe_both_designs_from_engine(composition_result, de_result):
    c = composition_result.methods
    assert (
        "a paired t-test on the per-donor differences of logit-transformed proportions when at "
        "least 2 donors have a sample in both groups, otherwise a Welch t-test on "
        "logit-transformed proportions"
    ) in c
    d = de_result.methods
    assert de_result.runs[0].details["design"] == "paired"
    assert (
        "~ pair + condition (blocked on the donor) when every pair has a pseudobulk in both "
        "groups, otherwise ~ condition"
    ) in d


def _data(**over):
    d = {
        "n_cells": 1800,
        "n_genes": 2500,
        "counts_source": ".X",
        "sha256": "0" * 64,
        "mito_prefix": "MT-",
        "paired_design": True,
    }
    d.update(over)
    return d


def _methods(spec, data):
    counts = {"one_at_a_time": 1, "random": 0}
    return methods_paragraph(spec, GridConfig(), Thresholds(), data, 2, counts, "none")


@pytest.mark.parametrize(
    "design,de_clause,comp_clause",
    [
        (
            "paired",
            "the baseline run used ~ pair + condition",
            "the baseline run used the paired test",
        ),
        (
            "unpaired",
            "the baseline run used ~ condition",
            "the baseline run used the unpaired test",
        ),
    ],
)
def test_methods_name_the_baseline_design(truth, design, de_clause, comp_clause):
    de = _spec("differential_expression", truth, comparison=DE, genes=["G0"])
    assert de_clause in _methods(de, _data(baseline_design=design))
    comp = _spec("composition_shift", truth, comparison=COMPOSITION)
    assert comp_clause in _methods(comp, _data(baseline_design=design))


def test_methods_omit_the_baseline_design_when_unknown(truth):
    de = _spec("differential_expression", truth, comparison=DE, genes=["G0"])
    assert "the baseline run used" not in _methods(de, _data())
    assert "the baseline run used" not in _methods(de, _data(baseline_design=None))
    # an unpaired dataset has one design, so there is nothing to tell apart
    text = _methods(de, _data(paired_design=False, baseline_design="unpaired"))
    assert "The design was ~ condition." in text
    assert "the baseline run used" not in text


# ------------------------------------------------------------ population panel


def test_population_report_shows_the_mixture_block(population_result):
    res = population_result
    dc = res.runs[0].match["doublet_cluster"]
    assert isinstance(dc, dict) and dc["best_pair"] and len(dc["best_pair"]) == 2
    html = render_report(res)
    assert "Cluster-level doublet check (baseline run)" in html
    assert "Best explaining pair" in html
    assert "Residual ratio, mixture vs best single cluster" in html
    assert "Library-size ratio, population vs mixture" in html
    assert "no cutoff on the residual ratio or the library-size ratio" in html
    assert "never changes the verdict" in html
    # the thresholds table says the two cutoffs are informational
    assert "doublet_libsize_ratio_min</span> are informational" in html


# -------------------------------------------------------------- time-of-day section


def _fake_result(truth, time_check, label_before=None):
    spec = _spec("composition_shift", truth, comparison=COMPOSITION)
    verdict = Verdict(
        label="FRAGILE",
        survival=0.5,
        n_runs=2,
        n_holds=1,
        baseline_holds=True,
        cell_level_survival=None,
        recovered_survival=None,
        pseudoreplication_flag=False,
        by_choice=[],
        breakers=[],
        oat_flips=[],
        summary="",
        label_before_time_check=label_before,
    )
    return ClaimResult(
        spec=spec,
        mode="run",
        verdict=verdict,
        runs=[],
        as_labelled=None,
        reference_source="none",
        data_summary={},
        grid=GridConfig(),
        thresholds=Thresholds(),
        provenance={},
        methods="",
        total_runtime_s=0.0,
        time_check=time_check,
    )


def test_time_section_renders_balance_and_adjustment_keys(truth):
    tc = {
        "mode": "A",
        "label": "POSSIBLE TIME CONFOUND",
        "reasons": ["collection time is associated with condition"],
        "disclaimer": "screen only",
        "balance": {
            "assessable": True,
            "unit_column": "donor",
            "n_samples_a": 6,
            "n_samples_b": 6,
            "n_units_a": 3,
            "n_units_b": 3,
            "circular_mean_hours_a": 9.0,
            "circular_mean_hours_b": 15.5,
            "circular_mean_difference_hours": -6.5,
            "mean_vector_distance": 1.2,
            "p_permutation": 0.1,
            "n_permutations": 20,
            "exact": True,
            "min_attainable_p": 0.05,
            "p_cos_ttest": 0.2,
            "p_sin_ttest": 0.3,
            "imbalanced": True,
            "imbalance_basis": "circular mean difference above 4 h",
        },
        "adjustment": {
            "design": "~ cos_t + sin_t + condition",
            "estimable": False,
            "reason": "time and condition are collinear across 3 donors",
            "held_unadjusted": True,
            "survives_adjustment": None,
            "shrunk": False,
            "flipped": False,
            "per_gene": [
                {
                    "gene": "G1",
                    "lfc_unadjusted": 1.2,
                    "padj_unadjusted": 0.01,
                    "lfc_adjusted": None,
                    "padj_adjusted": None,
                    "shrinkage": None,
                    "direction_flipped": False,
                    "survives": None,
                }
            ],
        },
        "panel": None,
        "rhythmic": None,
    }
    html = _time_html(_fake_result(truth, tc, label_before="ROBUST"))
    assert "POSSIBLE TIME CONFOUND" in html and "Mode A, collection time known" in html
    assert "Verdict changed from <strong>ROBUST</strong> to <strong>FRAGILE</strong>" in html
    assert "donor levels" in html and "3 vs 3" in html and "6 vs 6" in html
    assert "-6.50" in html
    assert "exact, 20 permutations; smallest attainable p 0.050" in html
    assert "yes (circular mean difference above 4 h)" in html
    assert "~ cos_t + sin_t + condition" in html
    assert "no: time and condition are collinear across 3 donors" in html
    assert "Held before adjustment</td><td class='num'>yes" in html
    assert html.count("not estimable") == 2  # the headline row and the per-gene row


def test_time_section_renders_panel_keys(truth):
    tc = {
        "mode": "B",
        "label": "NO TIME CONFOUND DETECTED",
        "reasons": [],
        "disclaimer": "d",
        "balance": None,
        "adjustment": None,
        "panel": {
            "assessable": True,
            "panel": ["PER1", "ARNTL"],
            "panel_citation": "cite",
            "genes_present": ["PER1", "ARNTL"],
            "genes_usable": ["PER1", "ARNTL"],
            "min_total_counts": 10,
            "unit_column": "donor",
            "n_samples_a": 4,
            "n_samples_b": 4,
            "index_genes": ["PER1", "ARNTL"],
            "index_difference_a_minus_b": 0.4,
            "index_cohens_d": 0.9,
            "index_p_permutation": 0.03,
            "index_n_permutations": 5000,
            "index_exact": False,
            "min_attainable_p": 0.0002,
            "index_welch_p": 0.04,
            "n_panel_genes_different": 1,
            "per_gene": [
                {
                    "gene": "PER1",
                    "arm": "ebox",
                    "mean_log2cpm_a": 5.0,
                    "mean_log2cpm_b": 4.0,
                    "diff_a_minus_b": 1.0,
                    "padj": 0.02,
                }
            ],
        },
        "rhythmic": {"available": False, "note": "no list"},
    }
    html = _time_html(_fake_result(truth, tc, label_before="FRAGILE"))
    assert "NO TIME CONFOUND DETECTED" in html and "Mode B, collection time unknown" in html
    # the label before the check is always recorded; an unchanged label is not flagged
    assert "Verdict changed" not in html
    assert "This check did not change the verdict." in html
    assert "donor levels" in html and "4 vs 4" in html
    assert "Genes in the index</td><td class='num'>PER1, ARNTL" in html
    assert "0.030 (5000 permutations; smallest attainable p 2.0e-04)" in html
    assert "Welch t on the index, p</td><td class='num'>0.040" in html
    assert "Per-gene panel values" in html and "ebox" in html
    assert "Day-night gene list: no list." in html


def test_time_section_tolerates_minimal_dicts(truth):
    tc = {
        "mode": "A",
        "label": "NOT ASSESSABLE",
        "reasons": ["no times"],
        "disclaimer": "d",
        "balance": {"assessable": True, "p_permutation": 0.5, "n_permutations": 10, "exact": True},
        "adjustment": {
            "effect_unadjusted": 0.3,
            "p_unadjusted": 0.2,
            "effect_adjusted": 0.1,
            "p_adjusted": 0.6,
            "shrinkage": 0.66,
            "survives_adjustment": False,
        },
        "panel": {"assessable": False, "reason": "too few genes"},
        "rhythmic": None,
    }
    html = _time_html(_fake_result(truth, tc))
    assert "NOT ASSESSABLE" in html and "exact, 10 permutations" in html
    assert "Composition refitted with time covariates" in html
    assert "Survives adjustment</td><td class='num'>no" in html
    assert "too few genes" in html
    assert "Not applicable" in _time_html(_fake_result(truth, None))


def test_time_section_in_engine_report(composition_result):
    html = render_report(composition_result)
    assert "Mode B, collection time unknown" in html
    assert "Clock-gene state across samples" in html
