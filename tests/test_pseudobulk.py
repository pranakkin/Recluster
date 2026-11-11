"""Pseudobulk aggregation must equal a plain per-sample sum and drop tiny samples."""

from __future__ import annotations

import numpy as np
import pandas as pd

from recluster.data import extract_counts
from recluster.evaluators.pseudobulk import pseudobulk_counts
from recluster.spec import make_sample_key


def test_pseudobulk_matches_manual_sum(sim):
    counts, _ = extract_counts(sim)
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cells = sim.obs.index[sim.obs["cell_type"] == "Mono"]
    pb, n_cells = pseudobulk_counts(counts, cells, sample, min_cells=1)

    dense = pd.DataFrame(counts[cells].X.toarray(), index=cells, columns=counts.var_names)
    manual = dense.groupby(sample.loc[cells].to_numpy()).sum()
    manual = manual.loc[pb.index]
    assert pb.shape == manual.shape
    assert np.array_equal(pb.to_numpy(), manual.to_numpy().astype(np.int64))
    assert pb.dtypes.unique().tolist() == [np.dtype("int64")]
    assert n_cells.sum() == len(cells)


def test_pseudobulk_drops_samples_below_min_cells(sim):
    counts, _ = extract_counts(sim)
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cells = sim.obs.index[sim.obs["cell_type"] == "Mono"]
    # keep only 3 cells from one sample
    small = sample.loc[cells] == sample.loc[cells].iloc[0]
    cells_small = cells[small.to_numpy()][:3].append(cells[~small.to_numpy()])
    pb, n_cells = pseudobulk_counts(counts, cells_small, sample, min_cells=5)
    assert sample.loc[cells].iloc[0] not in pb.index
    assert (n_cells >= 5).all()
    assert len(pb) == n_cells.shape[0]


def test_pseudobulk_total_counts_preserved(sim):
    counts, _ = extract_counts(sim)
    sample = make_sample_key(sim.obs, ["donor", "condition"])
    cells = sim.obs.index
    pb, _ = pseudobulk_counts(counts, cells, sample, min_cells=1)
    assert pb.to_numpy().sum() == int(counts.X.sum())
