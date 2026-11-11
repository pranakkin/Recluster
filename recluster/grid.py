"""Bounded perturbation design: baseline + one-at-a-time + seeded random combinations."""

from __future__ import annotations

import numpy as np

from recluster.config import GridConfig, PipelineParams


def build_grid(cfg: GridConfig, *, has_batch: bool) -> list[PipelineParams]:
    """Return the ordered list of runs. The baseline is always first.

    One-at-a-time runs change exactly one axis from baseline, so a failure there is
    attributable. Random combinations probe interactions. Duplicates are dropped.
    """
    base = cfg.baseline
    if not has_batch:
        base = base.model_copy(update={"batch_correction": False})
    axes: dict[str, list[object]] = dict(cfg.axes.model_dump())
    if not has_batch:
        axes.pop("batch_correction")

    runs: list[PipelineParams] = [base]
    for axis, values in axes.items():
        for v in values:
            if v == getattr(base, axis):
                continue
            p = base.model_copy(update={axis: v})
            if p not in runs:
                runs.append(p)

    rng = np.random.default_rng(cfg.random_seed)
    added, attempts = 0, 0
    while added < cfg.n_random and attempts < 50 * max(cfg.n_random, 1):
        attempts += 1
        update = {axis: values[int(rng.integers(len(values)))] for axis, values in axes.items()}
        p = base.model_copy(update=update)
        if p not in runs:
            runs.append(p)
            added += 1

    if cfg.max_runs is not None:
        runs = runs[: max(1, cfg.max_runs)]
    return runs


def describe_grid(runs: list[PipelineParams]) -> dict[str, int]:
    """Counts by run kind, for the report."""
    base = runs[0]
    n_oat = sum(1 for r in runs[1:] if len(r.diff(base)) == 1)
    return {
        "total": len(runs),
        "baseline": 1,
        "one_at_a_time": n_oat,
        "random": len(runs) - 1 - n_oat,
    }
