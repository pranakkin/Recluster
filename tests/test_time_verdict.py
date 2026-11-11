"""The time-of-day check folded into the verdict, and Mode B on planted clock-gene patterns."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from recluster.config import Thresholds
from recluster.results import Verdict
from recluster.spec import make_sample_key
from recluster.time_confound import (
    CLOCK_PANEL,
    INDEX_EXCLUDED,
    NO_CONFOUND,
    NOT_ASSESSABLE,
    POSSIBLE,
    decide,
    panel_check,
)
from recluster.verdict import FRAGILE, NOT_SUPPORTED, ROBUST, apply_time_check

# ------------------------------------------------------------------ apply_time_check

LABELS = (ROBUST, FRAGILE, NOT_SUPPORTED)
TIME_CHECKS = {
    "none": None,
    "a_removed": {"label": POSSIBLE, "mode": "A", "effect_removed": True, "headline": "x"},
    "a_shrunk": {"label": POSSIBLE, "mode": "A", "effect_shrunk": True},
    "b_possible": {"label": POSSIBLE, "mode": "B"},
    "no_confound": {"label": NO_CONFOUND},
    "not_assessable": {"label": NOT_ASSESSABLE},
}


def _verdict(label: str) -> Verdict:
    return Verdict(
        label=label,
        survival=1.0,
        n_runs=4,
        n_holds=4,
        baseline_holds=True,
        cell_level_survival=None,
        recovered_survival=None,
        pseudoreplication_flag=False,
        by_choice=[],
        breakers=[],
        oat_flips=[],
        summary="base.",
    )


def _expected_label(before: str, time_check: dict | None) -> str:
    """The documented rules: a removed effect sinks the claim, any other possible confound
    caps it at FRAGILE, and the two remaining labels leave it alone."""
    if time_check is None or time_check["label"] != POSSIBLE:
        return before
    if time_check.get("mode") == "A" and time_check.get("effect_removed"):
        return NOT_SUPPORTED
    return FRAGILE if before == ROBUST else before


@pytest.mark.parametrize("case", list(TIME_CHECKS))
@pytest.mark.parametrize("before", LABELS)
def test_apply_time_check_rules(before, case):
    time_check = TIME_CHECKS[case]
    v = apply_time_check(_verdict(before), time_check)
    assert v.label == _expected_label(before, time_check)
    if time_check is None:
        assert v.time_check_label is None and v.label_before_time_check is None
        assert v.summary == "base." and not v.changed_by_time_check
        return
    assert v.time_check_label == time_check["label"]
    assert v.label_before_time_check == before
    assert v.summary.startswith("base.") and "Time-of-day check" in v.summary
    assert v.changed_by_time_check == (v.label != before)


def test_apply_time_check_never_raises_a_label():
    for before in (FRAGILE, NOT_SUPPORTED):
        for case in ("a_shrunk", "b_possible", "no_confound", "not_assessable"):
            assert apply_time_check(_verdict(before), TIME_CHECKS[case]).label == before
    for before in LABELS:
        assert apply_time_check(_verdict(before), TIME_CHECKS["a_removed"]).label == NOT_SUPPORTED


def test_apply_time_check_summary_names_the_change():
    v = apply_time_check(_verdict(ROBUST), TIME_CHECKS["a_removed"])
    assert "changed from ROBUST to NOT SUPPORTED because x" in v.summary
    v = apply_time_check(_verdict(ROBUST), TIME_CHECKS["b_possible"])
    assert "changed from ROBUST to FRAGILE" in v.summary
    v = apply_time_check(_verdict(FRAGILE), TIME_CHECKS["b_possible"])
    assert "the verdict stays FRAGILE" in v.summary
    v = apply_time_check(_verdict(ROBUST), {"label": NO_CONFOUND, "headline": "balanced"})
    assert v.summary.endswith("no time confound detected (balanced).")


# ------------------------------------------------------------- Mode B on the sim

# Panel genes in the order they are assigned to the chosen sim genes (ranked by total counts),
# alternating arms so that both arms span the same expression range. CRY1 is in the panel but
# left out of the clock-state index.
GENES = ["PER1", "ARNTL", "PER2", "CLOCK", "PER3", "NPAS2", "CRY2", "NR1D1", "DBP", "TEF", "CRY1"]
ANTIPHASE = {g: (3.0 if CLOCK_PANEL[g] == "ebox" else 1 / 3) for g in GENES}
UNIFORM = {g: 3.0 for g in GENES}


def _with_clock_genes(sim: AnnData) -> AnnData:
    """The sim with 11 well-expressed background genes renamed to clock-panel genes.

    The background genes ranked 21-31 by total counts are used: high enough that per-sample
    log CPM is stable, below the heavy tail so that scaling them barely moves library sizes.
    """
    a = sim.copy()
    totals = pd.Series(np.asarray(a.X.sum(axis=0)).ravel(), index=a.var_names)
    background = totals[totals.index.str.startswith("G")].sort_values(ascending=False)
    chosen = list(background.index[20 : 20 + len(GENES)])
    names = a.var_names.to_series()
    names.loc[chosen] = GENES
    a.var_names = names.to_numpy()
    return a


def _scaled_in_stim(a: AnnData, factors: dict[str, float]) -> AnnData:
    """Multiply the counts of each gene in stim cells by its factor, rounded to integers."""
    b = a.copy()
    x = b.X.toarray()
    stim = (b.obs["condition"] == "stim").to_numpy()
    for gene, factor in factors.items():
        j = b.var_names.get_loc(gene)
        x[stim, j] = np.rint(x[stim, j] * factor)
    b.X = sp.csr_matrix(x.astype(np.float32))
    return b


def _panel(a: AnnData) -> dict:
    sample = make_sample_key(a.obs, ["donor", "condition"])
    group_by_sample = a.obs["condition"].astype(str).groupby(sample.to_numpy()).first()
    return panel_check(a, sample, group_by_sample, "stim", "ctrl", Thresholds())


@pytest.fixture(scope="module")
def clock_sim(sim):
    return _with_clock_genes(sim)


def test_panel_detects_an_antiphase_shift(clock_sim):
    panel = _panel(_scaled_in_stim(clock_sim, ANTIPHASE))
    assert panel["assessable"] is True
    assert panel["n_samples_a"] == 6 and panel["n_samples_b"] == 6
    assert set(panel["genes_usable"]) == set(GENES)
    assert "CRY1" in panel["genes_usable"] and "CRY1" not in panel["index_genes"]
    assert set(panel["index_genes"]) == set(GENES) - INDEX_EXCLUDED
    # 6 vs 6 units: the exact test has C(12, 6) = 924 assignments and a floor of 2/924
    assert panel["index_exact"] and panel["index_n_permutations"] == math.comb(12, 6) == 924
    assert panel["min_attainable_p"] == pytest.approx(2 / 924)
    assert panel["test_can_reject"] is True
    assert panel["index_p_permutation"] < 0.05
    assert panel["index_difference_a_minus_b"] > 0  # E-box arm up, activators down in stim
    assert panel["coherent_shift"] is True
    assert panel["n_panel_genes_different"] == len(GENES)
    tc = decide("B", None, None, panel, Thresholds())
    assert tc.label == POSSIBLE and tc.mode == "B"


def test_panel_uniform_shift_is_not_the_time_pattern(clock_sim):
    panel = _panel(_scaled_in_stim(clock_sim, UNIFORM))
    assert panel["assessable"] is True and panel["test_can_reject"] is True
    assert panel["any_gene_differs"] is True
    assert panel["n_panel_genes_different"] == len(GENES)
    assert panel["index_p_permutation"] >= 0.05
    assert panel["coherent_shift"] is False
    tc = decide("B", None, None, panel, Thresholds())
    assert tc.label == NO_CONFOUND
    assert any("not in the antiphase pattern" in r for r in tc.reasons)


def test_panel_three_vs_three_cannot_reject(clock_sim):
    a = _scaled_in_stim(clock_sim, ANTIPHASE)
    keep = a.obs["donor"].astype(str).isin(["donor0", "donor1", "donor2"]).to_numpy()
    panel = _panel(a[keep].copy())
    assert panel["assessable"] is True
    assert panel["n_samples_a"] == 3 and panel["n_samples_b"] == 3
    # C(6, 3) = 20 assignments, two of which tie at the maximum: the floor is 0.1
    assert panel["min_attainable_p"] == pytest.approx(2 / 20)
    assert panel["test_can_reject"] is False and panel["coherent_shift"] is False
    tc = decide("B", None, None, panel, Thresholds())
    assert tc.label == NOT_ASSESSABLE
    assert any("cannot reach" in r for r in tc.reasons)


def test_panel_check_is_repeatable(clock_sim):
    a = _scaled_in_stim(clock_sim, ANTIPHASE)
    first, second = _panel(a), _panel(a)
    assert first == second
    assert first["per_gene"] == second["per_gene"]
