"""Time-of-day confound check: balance, adjustment, and the not-assessable path."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from recluster.config import Thresholds
from recluster.hours import parse_hours
from recluster.spec import ClaimSpec, Comparison, make_sample_key
from recluster.time_confound import (
    NO_CONFOUND,
    NOT_ASSESSABLE,
    POSSIBLE,
    _balance_sentence,
    _permutation_p,
    adjusted_composition,
    adjusted_pseudobulk_de,
    assess_claim,
    balance_check,
    decide,
    exchangeable_units,
    panel_check,
    sample_times,
    summarise_adjustment,
)


def test_parse_hours_formats():
    s = pd.Series([9, 13.5, "08:30", "23:15:00", "2024-01-01 14:00", None, "nope"])
    out = parse_hours(s)
    assert out.iloc[0] == 9 and out.iloc[1] == 13.5
    assert math.isclose(out.iloc[2], 8.5) and math.isclose(out.iloc[3], 23.25)
    assert out.iloc[4] == 14 and math.isnan(out.iloc[5]) and math.isnan(out.iloc[6])


def _samples(n_per_group: int) -> pd.Series:
    return pd.Series(
        ["A"] * n_per_group + ["B"] * n_per_group,
        index=[f"s{i}" for i in range(2 * n_per_group)],
    )


def test_balance_balanced_vs_confounded():
    groups = _samples(6)
    balanced = pd.Series([8, 10, 12, 14, 16, 18, 8, 10, 12, 14, 16, 18], index=groups.index)
    confounded = pd.Series([7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19], index=groups.index)
    b = balance_check(balanced, groups, "A", "B")
    c = balance_check(confounded, groups, "A", "B")
    assert b["assessable"] and c["assessable"] and b["exact"] and c["exact"]
    assert b["p_permutation"] > 0.5 and not b["imbalanced"]
    assert c["p_permutation"] < 0.01 and c["imbalanced"]
    assert abs(c["circular_mean_difference_hours"] + 7) < 0.5  # A about 7 h earlier than B
    assert c["imbalance_basis"].startswith("permutation test and")
    assert b["test_can_reject"] and b["min_attainable_p"] == 2 / math.comb(12, 6)


def test_balance_with_few_samples_uses_the_hours_rule():
    # 3 vs 3: the exact test has 20 assignments and a floor of 0.1, so it can never reject
    groups = _samples(3)
    confounded = pd.Series([8, 8, 8, 16, 16, 16], index=groups.index, dtype=float)
    c = balance_check(confounded, groups, "A", "B")
    assert c["assessable"] and not c["test_can_reject"]
    assert c["min_attainable_p"] == pytest.approx(0.1) and c["p_permutation"] >= 0.1
    assert c["imbalanced"] and c["imbalance_basis"].startswith("mean difference")
    close = pd.Series([8, 9, 10, 9, 10, 11], index=groups.index, dtype=float)
    assert not balance_check(close, groups, "A", "B")["imbalanced"]  # 1 h apart


def test_balance_aligns_hours_and_groups_by_index():
    groups = _samples(6)
    hours = pd.Series(
        [7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19], index=groups.index, dtype=float
    )
    straight = balance_check(hours, groups, "A", "B")
    shuffled = balance_check(hours.iloc[::-1], groups, "A", "B")
    assert straight["circular_mean_difference_hours"] == shuffled["circular_mean_difference_hours"]


def test_sample_times_uses_the_circular_mean_across_midnight():
    cells = pd.Index([f"c{i}" for i in range(4)])
    times = pd.Series(["23:30", "00:30", "09:00", "09:00"], index=cells)
    sample = pd.Series(["night", "night", "day", "day"], index=cells)
    hours, warns = sample_times(times, sample)
    assert hours["night"] == pytest.approx(0.0) and hours["day"] == pytest.approx(9.0)
    assert warns and "more than one collection time" in warns[0]


def test_permutation_floor_and_monte_carlo_are_deterministic():
    rng = np.random.default_rng(0)
    values = rng.normal(size=24)
    is_a = np.array([True] * 12 + [False] * 12)

    def stat(v, m):
        return abs(float(v[m].mean() - v[~m].mean()))

    first = _permutation_p(values, is_a, stat)
    second = _permutation_p(values, is_a, stat)
    assert first == second and not first["exact"] and first["n_permutations"] == 5000
    small = _permutation_p(values[:6], np.array([True] * 3 + [False] * 3), stat)
    assert small["exact"] and small["min_attainable_p"] == pytest.approx(2 / 20)
    uneven = _permutation_p(values[:5], np.array([True] * 2 + [False] * 3), stat)
    assert uneven["min_attainable_p"] == pytest.approx(1 / 10)


def test_balance_needs_two_samples_per_group():
    groups = pd.Series(["A", "B", "B"], index=["s0", "s1", "s2"])
    hours = pd.Series([9.0, 10.0, 11.0], index=groups.index)
    assert balance_check(hours, groups, "A", "B")["assessable"] is False


def _pseudobulk_adata_with_time_effect(seed: int = 0) -> tuple[AnnData, pd.Series, pd.Series]:
    """Cells from 12 samples; one gene follows cos(time) and times are confounded with group."""
    rng = np.random.default_rng(seed)
    hours = {f"s{i}": h for i, h in enumerate([7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19])}
    group = {s: ("A" if i < 6 else "B") for i, s in enumerate(hours)}
    n_genes, cells_per = 300, 60
    base = rng.gamma(2.0, 2.0, size=n_genes)
    xs, obs_rows = [], []
    for s, h in hours.items():
        mult = np.ones(n_genes)
        mult[0] = 2 ** (1.5 * math.cos(2 * math.pi * (h - 7) / 24))  # peaks at 07:00
        lam = base * mult * rng.lognormal(0, 0.05, size=n_genes)
        xs.append(rng.poisson(lam[None, :] * rng.lognormal(0, 0.2, size=(cells_per, 1))))
        obs_rows += [{"sample": s, "group": group[s], "hour": h}] * cells_per
    x = np.vstack(xs).astype(np.float32)
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"c{i}" for i in range(len(obs))]
    a = AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))
    return a, a.obs["sample"].astype(str), a.obs["group"].astype(str)


def test_adjustment_removes_injected_time_effect():
    a, sample, group = _pseudobulk_adata_with_time_effect()
    thr = Thresholds()
    hours = a.obs.groupby("sample")["hour"].first().astype(float)
    cmp = Comparison(group_a="A", group_b="B", direction="up_in_a")
    df, info = adjusted_pseudobulk_de(
        a, a.obs_names, sample, group, hours, cmp, ["G0"], thr, n_cpus=2
    )
    assert info["estimable"] and info["design"] == "~condition"
    summary = summarise_adjustment(df, "all", thr, direction="up_in_a", info=info)
    g0 = next(r for r in summary["per_gene"] if r["gene"] == "G0")
    assert g0["held_unadjusted"], g0  # the confounded design calls a 'group effect'
    assert g0["survives"] is False, g0  # which the time covariates remove
    assert g0["shrinkage"] is not None and g0["shrinkage"] > 0.5
    tc = decide(
        "A",
        {"assessable": True, "p_permutation": 0.001, "circular_mean_difference_hours": -7.0},
        summary,
        None,
        thr,
    )
    assert tc.label == POSSIBLE and tc.effect_removed
    # the claimed direction matters: the same fit read as 'down_in_a' never held
    down = summarise_adjustment(df, "all", thr, direction="down_in_a", info=info)
    assert down["held_unadjusted"] is False


def test_refit_is_skipped_when_time_is_collinear_with_condition():
    a, sample, group = _pseudobulk_adata_with_time_effect()
    thr = Thresholds()
    cmp = Comparison(group_a="A", group_b="B", direction="either")
    two_valued = pd.Series(
        [9.0 if g == "A" else 15.5 for g in a.obs.groupby("sample")["group"].first()],
        index=a.obs.groupby("sample")["group"].first().index,
    )
    df, info = adjusted_pseudobulk_de(
        a, a.obs_names, sample, group, two_valued, cmp, ["G0"], thr, n_cpus=2
    )
    assert info["estimable"] is False and "collinear" in info["reason"]
    assert df["p_adjusted"].isna().all()
    summary = summarise_adjustment(df, "all", thr, info=info)
    assert summary["survives_adjustment"] is None
    tc = decide(
        "A",
        {"assessable": True, "p_permutation": 0.001, "circular_mean_difference_hours": -6.5},
        summary,
        None,
        thr,
    )
    assert tc.label == POSSIBLE and not tc.effect_removed
    assert any("cannot be separated" in r for r in tc.reasons)
    # too few pseudobulks for four model terms
    few = a[a.obs["sample"].isin(["s0", "s1", "s6", "s7"]).to_numpy()]
    df2, info2 = adjusted_pseudobulk_de(
        few,
        few.obs_names,
        sample.loc[few.obs_names],
        group.loc[few.obs_names],
        a.obs.groupby("sample")["hour"].first().astype(float),
        cmp,
        ["G0"],
        thr,
        n_cpus=2,
    )
    assert info2["estimable"] is False and "too few" in info2["reason"]


def test_composition_adjustment_filters_groups_and_reports_rank_deficiency():
    rng = np.random.default_rng(1)
    idx = [f"s{i}" for i in range(9)]
    table = pd.DataFrame(
        {
            "condition": ["A"] * 3 + ["B"] * 3 + ["C"] * 3,
            "logit": np.r_[
                rng.normal(1.0, 0.1, 3), rng.normal(0.0, 0.1, 3), rng.normal(-3, 0.1, 3)
            ],
        },
        index=idx,
    )
    hours = pd.Series([8, 9, 10, 14, 15, 16, 20, 21, 22], index=idx, dtype=float)
    adj = adjusted_composition(table, hours, "A", "B")
    assert adj["n_samples_a"] == 3 and adj["n_samples_b"] == 3 and adj["estimable"]
    assert adj["effect_unadjusted"] == pytest.approx(1.0, abs=0.3)
    two_valued = pd.Series([9] * 3 + [15] * 3 + [20] * 3, index=idx, dtype=float)
    bad = adjusted_composition(table, two_valued, "A", "B")
    assert bad["estimable"] is False and math.isnan(bad["p_adjusted"])


def test_exchangeable_units_collapse_nested_samples():
    donors = [f"d{i}" for i in range(8)]
    rows = []
    for i, d in enumerate(donors):
        for lab in ("ctrl", "stim"):
            rows += [{"donor": d, "label": lab, "group": "A" if i < 4 else "B"}] * 5
    obs = pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))])
    sample = make_sample_key(obs, ["donor", "label"])
    cond = obs["group"]
    unit, col, notes = exchangeable_units(
        obs, obs.index, sample, cond, "A", "B", ["donor", "label"], "group"
    )
    assert col == "donor" and unit.nunique() == 8 and notes and "nested" in notes[0]
    # condition varies within the donor: the sample stays the unit
    unit2, col2, _ = exchangeable_units(
        obs, obs.index, sample, obs["label"], "stim", "ctrl", ["donor", "label"], "label"
    )
    assert col2 is None and unit2.nunique() == 16
    # a donor-level time that differs between a donor's samples blocks the collapse
    hours = pd.Series(
        [9.0 if s.endswith("ctrl") else 15.0 for s in sample.unique()], index=sample.unique()
    )
    _, col3, _ = exchangeable_units(
        obs, obs.index, sample, cond, "A", "B", ["donor", "label"], "group", hours_by_sample=hours
    )
    assert col3 is None


def test_decide_rules_cover_each_branch():
    thr = Thresholds()
    bal_ok = {"assessable": True, "p_permutation": 0.6, "circular_mean_difference_hours": 0.5}
    assert decide("A", bal_ok, None, None, thr).label == NO_CONFOUND
    bal_bad = {"assessable": True, "p_permutation": 0.01, "circular_mean_difference_hours": -6.0}
    survives = {
        "held_unadjusted": True,
        "survives_adjustment": True,
        "shrunk": False,
        "flipped": False,
    }
    assert decide("A", bal_bad, survives, None, thr).label == NO_CONFOUND
    shrunk = dict(survives, shrunk=True)
    tc = decide("A", bal_bad, shrunk, None, thr)
    assert tc.label == POSSIBLE and tc.effect_shrunk and not tc.effect_removed
    # imbalanced but the claim never held: possible confound, nothing removed
    not_held = dict(survives, held_unadjusted=False)
    tc = decide("A", bal_bad, not_held, None, thr)
    assert tc.label == POSSIBLE and not tc.effect_removed
    # the baseline run's own outcome wins over the refit's
    tc = decide("A", bal_bad, not_held, None, thr, held_unadjusted=True)
    assert tc.label == NO_CONFOUND
    # not significant after adjustment with a shrunk estimate: removed
    removed = dict(survives, survives_adjustment=False, shrunk=True)
    tc = decide("A", bal_bad, removed, None, thr)
    assert tc.label == POSSIBLE and tc.effect_removed
    # not significant but the estimate barely moved: cannot separate, not removed
    df_loss = dict(survives, survives_adjustment=False)
    tc = decide("A", bal_bad, df_loss, None, thr)
    assert tc.label == POSSIBLE and not tc.effect_removed
    # undefined adjusted estimate
    undefined = dict(survives, survives_adjustment=None)
    assert decide("A", bal_bad, undefined, None, thr).label == POSSIBLE
    # not estimable design
    tc = decide("A", bal_bad, {"estimable": False, "reason": "collinear"}, None, thr)
    assert tc.label == POSSIBLE and not tc.effect_removed
    # imbalanced by the hours rule alone with no adjustment possible
    bal_hours = {"assessable": True, "p_permutation": 0.1, "circular_mean_difference_hours": 8.0}
    assert decide("A", bal_hours, None, None, thr).label == POSSIBLE
    assert (
        decide("B", None, None, {"assessable": False, "reason": "x"}, thr).label == NOT_ASSESSABLE
    )
    underpowered = {
        "assessable": True,
        "test_can_reject": False,
        "min_attainable_p": 0.1,
        "n_samples_a": 3,
        "n_samples_b": 3,
        "index_difference_a_minus_b": 1.2,
        "index_welch_p": 0.01,
        "genes_usable": list("abcdef"),
        "n_panel_genes_different": 2,
        "genes_different": ["a", "b"],
        "coherent_shift": False,
        "any_gene_differs": True,
        "index_p_permutation": 0.1,
    }
    tc = decide("B", None, None, underpowered, thr)
    assert tc.label == NOT_ASSESSABLE and "cannot reach" in tc.reasons[0]
    coherent = dict(
        underpowered, test_can_reject=True, coherent_shift=True, index_p_permutation=0.01
    )
    assert decide("B", None, None, coherent, thr).label == POSSIBLE
    uniform = dict(
        underpowered, test_can_reject=True, coherent_shift=False, index_p_permutation=0.4
    )
    tc = decide("B", None, None, uniform, thr)
    assert tc.label == NO_CONFOUND and "not in the antiphase pattern" in tc.reasons[0]
    assert decide("none", None, None, None, thr).label == NOT_ASSESSABLE


def test_mode_b_not_assessable_without_clock_genes(sim):
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    group_by_sample = sim.obs["condition"].astype(str).groupby(sample.to_numpy()).first()
    panel = panel_check(sim, sample, group_by_sample, "stim", "ctrl", Thresholds())
    assert panel["assessable"] is False
    assert panel["genes_usable"] == []
    spec = ClaimSpec.model_validate(
        {
            "claim_type": "composition_shift",
            "name": "x",
            "columns": {"condition": "condition", "sample": ["donor", "condition"]},
            "target": {"markers": ["Mono_marker_0"]},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "different"},
        }
    )
    tc = assess_claim(
        spec,
        sim,
        sample,
        sim.obs["condition"].astype(str),
        sim.obs,
        None,
        sim.obs_names,
        Thresholds(),
        n_cpus=2,
    )
    assert tc.mode == "B" and tc.label == NOT_ASSESSABLE
    assert "panel genes" in tc.reasons[0]


def test_spec_validates_collection_time_column(sim):
    a = sim.copy()
    a.obs["collected"] = np.where(a.obs["condition"] == "ctrl", "09:00", "15:30")
    spec = ClaimSpec.model_validate(
        {
            "claim_type": "composition_shift",
            "name": "x",
            "columns": {
                "condition": "condition",
                "sample": ["donor", "condition"],
                "collection_time": "collected",
            },
            "target": {"markers": ["Mono_marker_0"]},
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "different"},
        }
    )
    spec.check_against(a)
    a.obs["collected"] = "sometime"
    with pytest.raises(Exception, match="collection_time"):
        spec.check_against(a)


def test_two_clock_times_balanced_across_groups_keep_one_time_term():
    a, sample, group = _pseudobulk_adata_with_time_effect()
    thr = Thresholds()
    cmp = Comparison(group_a="A", group_b="B", direction="either")
    per_sample = a.obs.groupby("sample")["group"].first()
    # 09:00 and 13:00 both occur in both groups: with two clock times sin is an affine function
    # of cos, so one time term is kept and the condition stays estimable
    hours = pd.Series(
        [9.0 if i % 2 == 0 else 13.0 for i in range(len(per_sample))], index=per_sample.index
    )
    df, info = adjusted_pseudobulk_de(
        a, a.obs_names, sample, group, hours, cmp, ["G0"], thr, n_cpus=2
    )
    assert info["estimable"] and info["time_terms"] == ["cos_t"]
    assert df["p_adjusted"].notna().all()


def test_balance_sentence_names_the_groups() -> None:
    s = _balance_sentence(
        {
            "group_a": "B",
            "group_b": "A",
            "circular_mean_difference_hours": 5.8,
            "p_permutation": 0.03,
        }
    )
    assert "+5.8 h (B minus A)" in s
