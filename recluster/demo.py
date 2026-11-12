"""A short demonstration on a seeded subset of Kang et al. 2018.

Two claims run on the full default grid: one true (ISG15 rises in CD14+ monocytes after
IFN-beta) and one false by construction (ISG15 differs between two fake groups of control
donors). The subset keeps all 16 donor x condition samples and draws the same number of cells
from each, so the replicate structure is intact while each run is fast.
"""

from __future__ import annotations

import os
import time
import webbrowser
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from anndata import AnnData

from recluster.bench.kang import (
    COLS,
    MONO,
    download_kang,
    fake_split_within_control,
    load_kang,
)
from recluster.config import GridConfig
from recluster.spec import ClaimSpec

CELLS_PER_SAMPLE = 200
SEED = 0


@dataclass
class DemoResult:
    name: str
    expected: str
    verdict: str
    n_holds: int
    n_runs: int
    report: Path
    seconds: float


def demo_subset(adata: AnnData, *, per_sample: int = CELLS_PER_SAMPLE, seed: int = SEED) -> AnnData:
    """The same number of random cells from each donor x condition sample (all if fewer)."""
    rng = np.random.default_rng(seed)
    sample = (adata.obs["replicate"].astype(str) + "|" + adata.obs["label"].astype(str)).to_numpy()
    keep: list[int] = []
    for s in sorted(set(sample)):
        idx = np.flatnonzero(sample == s)
        take = idx if len(idx) <= per_sample else rng.choice(idx, size=per_sample, replace=False)
        keep.extend(int(i) for i in take)
    return adata[np.sort(keep)].copy()


def demo_claims() -> tuple[ClaimSpec, ClaimSpec]:
    true_claim = ClaimSpec.model_validate(
        {
            "claim_type": "differential_expression",
            "name": "ISG15 up in CD14+ monocytes after IFN-beta (Kang 2018 subset)",
            "columns": COLS,
            "target": MONO,
            "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "up_in_a"},
            "genes": ["ISG15"],
        }
    )
    false_claim = ClaimSpec.model_validate(
        {
            "claim_type": "differential_expression",
            "name": "FALSE: ISG15 differs between fake donor groups (Kang subset, control cells)",
            "columns": {
                "condition": "fake_group",
                "sample": ["replicate", "fake_group"],
                "labels": "cell_type",
            },
            "target": MONO,
            "comparison": {"group_a": "A", "group_b": "B", "direction": "either"},
            "genes": ["ISG15"],
        }
    )
    return true_claim, false_claim


def default_jobs() -> int:
    """Worker processes for the demo: all but one core, at most 8."""
    return max(1, min(8, (os.cpu_count() or 2) - 1))


def run_demo(
    *,
    data_dir: Path = Path("data"),
    out_dir: Path = Path("results/demo"),
    jobs: int | None = None,
    grid: GridConfig | None = None,
    open_reports: bool = True,
    say=print,  # noqa: ANN001
) -> list[DemoResult]:
    """Download Kang 2018 if needed, build the subset and run both claims on the default grid.

    The claims run one after the other, each spreading its grid runs over ``jobs`` processes;
    running the two claims side by side was measured to be slower (nested process pools).
    """
    from recluster.engine import run_claim, slugify

    jobs = jobs or default_jobs()
    grid = grid or GridConfig()
    out_dir.mkdir(parents=True, exist_ok=True)
    say(f"Kang et al. 2018: {data_dir / 'kang_2018.h5ad'} (38 MB download on first use)")
    full = load_kang(download_kang(data_dir / "kang_2018.h5ad"))
    subset = demo_subset(full)
    split = fake_split_within_control(subset, seed=SEED)
    subset_path = out_dir / "kang_demo_subset.h5ad"
    split_path = out_dir / "kang_demo_fake_split.h5ad"
    subset.write_h5ad(subset_path)
    split.write_h5ad(split_path)
    say(
        f"subset: {subset.n_obs:,} cells ({CELLS_PER_SAMPLE} per donor x condition sample, "
        f"seed {SEED}); fake groups A = {', '.join(split.uns['fake_split']['A'])}"
    )
    true_claim, false_claim = demo_claims()
    results = []
    for data_path, expected, spec in (
        (subset_path, "ROBUST", true_claim),
        (split_path, "NOT SUPPORTED", false_claim),
    ):
        say(f"running on the default grid with {jobs} worker(s): {spec.name}")
        t0 = time.perf_counter()
        res = run_claim(data_path, spec, out_dir=out_dir, grid=grid, progress=None, jobs=jobs)
        v = res.verdict
        r = DemoResult(
            name=spec.name,
            expected=expected,
            verdict=v.label,
            n_holds=v.n_holds,
            n_runs=v.n_runs,
            report=out_dir / slugify(spec.name) / "report.html",
            seconds=time.perf_counter() - t0,
        )
        results.append(r)
        say(
            f"  {r.verdict} (expected {r.expected}): held in {r.n_holds} of {r.n_runs} runs, "
            f"{r.seconds:.0f} s\n  {r.report}"
        )
        if open_reports:
            webbrowser.open(r.report.resolve().as_uri())
    return results
