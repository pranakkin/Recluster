"""Sum raw counts per biological replicate within a set of cells."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


def pseudobulk_counts(
    counts: AnnData,
    cells: pd.Index,
    sample: pd.Series,
    *,
    min_cells: int,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (samples x genes integer counts, cells per sample).

    Samples contributing fewer than ``min_cells`` cells are dropped: a pseudobulk made of two
    cells is not a replicate.
    """
    idx = counts.obs_names.get_indexer(cells)
    idx = idx[idx >= 0]
    groups = sample.loc[counts.obs_names[idx]].astype(str).to_numpy()
    uniq, inverse = np.unique(groups, return_inverse=True)
    indicator = sp.csr_matrix(
        (np.ones(len(idx), dtype=np.float64), (inverse, np.arange(len(idx)))),
        shape=(len(uniq), len(idx)),
    )
    summed = indicator @ counts.X[idx]
    summed = summed.toarray() if sp.issparse(summed) else np.asarray(summed)
    df = pd.DataFrame(np.rint(summed).astype(np.int64), index=uniq, columns=counts.var_names)
    n_cells = pd.Series(np.bincount(inverse, minlength=len(uniq)), index=uniq, name="n_cells")
    keep = n_cells >= min_cells
    return df.loc[keep], n_cells.loc[keep]
