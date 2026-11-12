"""The demo's subset and claims; the download and full runs are exercised by `recluster demo`."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from typer.testing import CliRunner

from recluster.cli import app
from recluster.demo import default_jobs, demo_claims, demo_subset


def _kang_shaped(n_per_sample: int = 60) -> AnnData:
    rng = np.random.default_rng(0)
    rows = []
    for d in range(4):
        for label in ("ctrl", "stim"):
            n = n_per_sample if d else n_per_sample // 3  # one donor has fewer cells
            rows += [{"replicate": f"d{d}", "label": label, "cell_type": "CD14+ Monocytes"}] * n
    obs = pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))])
    x = sp.csr_matrix(rng.poisson(1.0, size=(len(obs), 20)).astype(np.float32))
    return AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(20)]))


def test_subset_draws_the_same_number_from_each_sample():
    a = _kang_shaped()
    sub = demo_subset(a, per_sample=30, seed=0)
    per = sub.obs.groupby(["replicate", "label"], observed=True).size()
    assert per.loc[("d1", "ctrl")] == 30 and per.loc[("d3", "stim")] == 30
    assert per.loc[("d0", "ctrl")] == 20  # fewer than requested: all kept
    again = demo_subset(a, per_sample=30, seed=0)
    assert list(again.obs_names) == list(sub.obs_names)
    assert list(demo_subset(a, per_sample=30, seed=1).obs_names) != list(sub.obs_names)


def test_demo_claims_are_one_true_and_one_false_by_construction():
    true_claim, false_claim = demo_claims()
    assert true_claim.comparison is not None and false_claim.comparison is not None
    assert true_claim.comparison.group_a == "stim" and true_claim.genes == ["ISG15"]
    assert false_claim.columns.condition == "fake_group"
    assert false_claim.columns.sample_key == ["replicate", "fake_group"]
    assert false_claim.name.startswith("FALSE")


def test_demo_command_is_listed_and_jobs_default_is_bounded():
    out = CliRunner().invoke(app, ["demo", "--help"]).output
    assert "--no-open" in out and "--jobs" in out
    assert 1 <= default_jobs() <= 8
