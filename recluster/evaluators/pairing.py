"""Detect a paired design: the same biological unit observed in both comparison groups."""

from __future__ import annotations

import pandas as pd


def pair_key(
    sample: pd.Series,
    condition: pd.Series,
    cells: pd.Index,
    group_a: str,
    group_b: str,
    sample_columns: list[str] | None,
    condition_column: str | None,
    obs: pd.DataFrame,
) -> pd.Series | None:
    """Return a per-cell pair id when every pair level appears in both groups, else None.

    The pair id is the replicate key without the condition column (for Kang 2018:
    ``[replicate, label]`` minus ``label`` gives the donor). Only the cells of the two compared
    groups count: a donor observed under some third condition alone does not make the design
    unpaired, and the returned series is indexed by the cells of the two groups only. Pairing
    is declared when every pair level has at least one sample in group A and one in group B
    and there are at least two levels; a partially paired design falls back to the unpaired
    tests.
    """
    if not sample_columns or condition_column is None or condition_column not in sample_columns:
        return None
    pair_cols = [c for c in sample_columns if c != condition_column]
    if not pair_cols:
        return None
    cond = condition.loc[cells].astype(str)
    in_groups = cond.isin([group_a, group_b]).to_numpy()
    cells, cond = cells[in_groups], cond[in_groups]
    if len(cells) == 0:
        return None
    sub = obs.loc[cells, pair_cols].astype(str)
    key = sub.iloc[:, 0]
    for c in pair_cols[1:]:
        key = key + "|" + sub[c]
    table = pd.crosstab(key.to_numpy(), cond.to_numpy())
    if group_a not in table.columns or group_b not in table.columns:
        return None
    if ((table[group_a] > 0) & (table[group_b] > 0)).all() and len(table) >= 2:
        return pd.Series(key.to_numpy(), index=cells, name="pair")
    return None
