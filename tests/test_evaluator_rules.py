"""Evaluator preconditions: pairing ignores conditions outside the comparison, composition needs
the population in enough samples of one group only, DE still needs it in both groups."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from recluster.config import Thresholds
from recluster.data import extract_counts
from recluster.evaluators.composition import evaluate_composition, group_presence
from recluster.evaluators.de import evaluate_de
from recluster.evaluators.pairing import pair_key
from recluster.matching import lognorm_matrix
from recluster.spec import Comparison, make_sample_key

SAMPLE_COLS = ["donor", "condition"]
STIM_VS_CTRL = {"group_a": "stim", "group_b": "ctrl"}


def _pair(obs: pd.DataFrame, cells: pd.Index) -> pd.Series | None:
    sample = make_sample_key(obs, SAMPLE_COLS)
    cond = obs["condition"].astype(str)
    return pair_key(sample, cond, cells, "stim", "ctrl", SAMPLE_COLS, "condition", obs)


def _mono(sim, condition: str | None = None) -> pd.Index:
    mask = (sim.obs["cell_type"] == "Mono").to_numpy()
    if condition is not None:
        mask &= (sim.obs["condition"] == condition).to_numpy()
    return sim.obs_names[mask]


def _relabel(sim, mask: np.ndarray) -> pd.DataFrame:
    """obs with the condition of the masked cells set to a third value."""
    obs = sim.obs.copy()
    obs["condition"] = np.where(mask, "other", obs["condition"].astype(str).to_numpy())
    return obs


@pytest.fixture(scope="module")
def third(sim) -> pd.DataFrame:
    """Every cell of donor0 observed under a third condition."""
    return _relabel(sim, (sim.obs["donor"].astype(str) == "donor0").to_numpy())


@pytest.fixture(scope="module")
def ctx(sim):
    counts, _ = extract_counts(sim)
    genes = list(dict(sim.uns["truth"])["de_genes"][:2])
    return {"counts": counts, "genes": genes, "lognorm": lognorm_matrix(counts, genes)}


# ------------------------------------------------------------------ pairing


def test_pairing_ignores_a_donor_seen_only_under_a_third_condition(sim, third):
    full = _pair(sim.obs, sim.obs_names)
    assert full is not None and full.nunique() == 6 and len(full) == sim.n_obs
    pair = _pair(third, sim.obs_names)
    assert pair is not None and pair.nunique() == 5
    assert "donor0" not in set(pair)
    # the pair id exists for the cells of the two compared groups only
    in_groups = third["condition"].isin(["ctrl", "stim"])
    assert len(pair) == int(in_groups.sum())
    assert set(third.loc[pair.index, "condition"]) == {"ctrl", "stim"}


def test_pairing_still_needs_every_remaining_donor_in_both_groups(sim):
    donor = sim.obs["donor"].astype(str).to_numpy()
    cond = sim.obs["condition"].astype(str).to_numpy()
    obs = _relabel(sim, (donor == "donor0") & (cond == "stim"))
    assert _pair(obs, sim.obs_names) is None  # donor0 keeps its ctrl cells but has no stim cells


def test_evaluators_accept_a_pair_key_restricted_to_the_compared_groups(sim, third, ctx):
    sample = make_sample_key(third, SAMPLE_COLS)
    cond = third["condition"].astype(str)
    pair = _pair(third, sim.obs_names)
    thr = Thresholds()
    # the run and the population include donor0's cells, which have no pair id
    res = evaluate_composition(
        _mono(sim),
        sim.obs_names,
        sample,
        cond,
        Comparison(direction="higher_in_a", **STIM_VS_CTRL),
        thr,
        pair=pair,
    )
    assert res["design"] == "paired" and res["holds"] and res["evaluable"]
    assert res["n_samples_a"] == 5 and res["n_samples_b"] == 5
    other = [r for r in res["per_sample"] if r["condition"] == "other"]
    assert len(other) == 1 and pd.isna(other[0]["pair"])

    de = evaluate_de(
        cells=_mono(sim),
        counts=ctx["counts"],
        lognorm=ctx["lognorm"],
        sample=sample,
        condition=cond,
        comparison=Comparison(direction="up_in_a", **STIM_VS_CTRL),
        genes=ctx["genes"],
        require="all",
        thr=thr,
        n_cpus=2,
        pair=pair,
    )
    assert de["design"] == "paired" and de["design_formula"] == "~pair + condition"
    assert de["n_samples_a"] == 5 and de["n_samples_b"] == 5
    assert de["holds"] and de["evaluable"]


# -------------------------------------------------------------- composition


def test_group_presence_needs_enough_samples_in_both_groups():
    thr = Thresholds()
    cmp = Comparison(direction="different", **STIM_VS_CTRL)
    table = pd.DataFrame(
        {"n": [100] * 4, "k": [30, 25, 10, 7], "condition": ["stim", "stim", "ctrl", "ctrl"]},
        index=["a1", "a2", "b1", "b2"],
    )
    assert group_presence(table, cmp, thr) == (True, [])
    table["k"] = [30, 25, 0, 2]  # the population is absent from the ctrl samples
    ok, reasons = group_presence(table, cmp, thr)
    assert not ok and len(reasons) == 1 and "'ctrl'" in reasons[0]
    assert "0 sample(s)" in reasons[0] and "cannot test" in reasons[0]
    table["k"] = [30, 2, 0, 7]  # one sample per group holds it
    ok, reasons = group_presence(table, cmp, thr)
    assert not ok and len(reasons) == 2


def test_composition_population_absent_from_one_group_is_not_testable(sim):
    # absence from one condition cannot be told apart from the clustering having lost that
    # condition's version of the population, so the run is reported as not testable
    sample = make_sample_key(sim.obs, SAMPLE_COLS)
    cond = sim.obs["condition"].astype(str)
    cmp = Comparison(direction="lower_in_a", **STIM_VS_CTRL)
    thr = Thresholds()
    res = evaluate_composition(_mono(sim, "ctrl"), sim.obs_names, sample, cond, cmp, thr)
    assert res["evaluable"] is False and res["holds"] is False
    assert len(res["reasons"]) == 1 and "'stim'" in res["reasons"][0]
    assert res["holds_cell_level"]  # the cell-level test would have called it
    assert res["mean_proportion_a"] == 0.0 and res["n_samples_a"] == 6
    assert all(r["k"] == 0 for r in res["per_sample"] if r["condition"] == "stim")
    paired = evaluate_composition(
        _mono(sim, "ctrl"),
        sim.obs_names,
        sample,
        cond,
        cmp,
        thr,
        pair=_pair(sim.obs, sim.obs_names),
    )
    assert paired["design"] == "paired" and paired["evaluable"] is False
    assert paired["design_note"] is None


def test_composition_not_testable_with_the_population_in_one_sample(sim):
    sample = make_sample_key(sim.obs, SAMPLE_COLS)
    cond = sim.obs["condition"].astype(str)
    cmp = Comparison(direction="lower_in_a", **STIM_VS_CTRL)
    one_sample = sim.obs_names[
        ((sim.obs["cell_type"] == "Mono") & (sim.obs["sample"] == "donor0|ctrl")).to_numpy()
    ]
    res = evaluate_composition(one_sample, sim.obs_names, sample, cond, cmp, Thresholds())
    assert res["evaluable"] is False and res["holds"] is False
    assert len(res["reasons"]) == 2 and all("cannot test" in r for r in res["reasons"])
    assert all(r in res["warnings"] for r in res["reasons"])


def test_composition_not_testable_with_one_replicate_in_a_group(sim):
    sample = make_sample_key(sim.obs, SAMPLE_COLS)
    cond = sim.obs["condition"].astype(str)
    cmp = Comparison(direction="higher_in_a", **STIM_VS_CTRL)
    # the population is present in every ctrl sample, but the run holds a single stim sample
    run_cells = sim.obs_names[
        ((sim.obs["condition"] == "ctrl") | (sim.obs["sample"] == "donor0|stim")).to_numpy()
    ]
    pop = _mono(sim).intersection(run_cells)
    res = evaluate_composition(pop, run_cells, sample, cond, cmp, Thresholds())
    assert res["n_samples_a"] == 1 and res["n_samples_b"] == 6
    assert res["evaluable"] is False and res["holds"] is False
    assert res["p_logit_ttest"] is None
    assert any("fewer than 2 replicates per group (1 vs 6)" in r for r in res["reasons"])


# ------------------------------------------------------------------------ DE


def _de(ctx, cells, sample, condition, pair=None):
    return evaluate_de(
        cells=cells,
        counts=ctx["counts"],
        lognorm=ctx["lognorm"],
        sample=sample,
        condition=condition,
        comparison=Comparison(direction="up_in_a", **STIM_VS_CTRL),
        genes=ctx["genes"],
        require="all",
        thr=Thresholds(),
        n_cpus=2,
        pair=pair,
    )


def test_de_still_needs_the_population_in_both_groups(sim, ctx):
    sample = make_sample_key(sim.obs, SAMPLE_COLS)
    cond = sim.obs["condition"].astype(str)
    res = _de(ctx, _mono(sim, "ctrl"), sample, cond, pair=_pair(sim.obs, sim.obs_names))
    assert res["evaluable"] is False and res["holds"] is False
    assert res["n_samples_a"] == 0 and res["n_samples_b"] == 6
    assert len(res["reasons"]) == 1 and "cannot test" in res["reasons"][0]
    assert "0 vs 6 replicate(s)" in res["reasons"][0]
    # nothing was fitted: no formula, no fallback
    assert res["design"] == "unpaired" and res["design_formula"] is None
    assert res["design_note"] is None


def test_de_reports_the_design_formula_and_the_fallback(sim, ctx):
    sample = make_sample_key(sim.obs, SAMPLE_COLS)
    cond = sim.obs["condition"].astype(str)
    pair = _pair(sim.obs, sim.obs_names)
    res = _de(ctx, _mono(sim), sample, cond, pair=pair)
    assert res["design"] == "paired" and res["design_formula"] == "~pair + condition"
    assert res["design_note"] is None and res["holds"]
    # without donor0's stim cells donor0 has a pseudobulk in ctrl only: the design falls back
    mono = _mono(sim)
    cells = mono[~(sim.obs.loc[mono, "sample"] == "donor0|stim").to_numpy()]
    res = _de(ctx, cells, sample, cond, pair=pair)
    assert res["design"] == "unpaired" and res["design_formula"] == "~condition"
    assert res["n_samples_a"] == 5 and res["n_samples_b"] == 6
    assert (
        res["design_note"] == "not every pair has a pseudobulk in both groups; unpaired design used"
    )
    assert res["design_note"] in res["warnings"]
