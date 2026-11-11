"""Typed containers for per-run and per-claim results."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from recluster.config import GridConfig, Thresholds
from recluster.spec import ClaimSpec


@dataclass
class RunEvaluation:
    """One grid run: pipeline bookkeeping + where the population landed + the test outcome."""

    index: int
    label: str
    params: dict[str, Any]
    is_baseline: bool
    is_one_at_a_time: bool
    n_cells: int
    n_removed_qc: int
    n_removed_doublet: int
    n_clusters: int
    runtime_s: float
    warnings: list[str]
    match: dict[str, Any]
    holds: bool
    holds_cell_level: bool | None
    replicate_ok: bool | None
    recovered: bool | None
    key_stats: dict[str, float | None]
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Verdict:
    label: str  # ROBUST | FRAGILE | NOT SUPPORTED
    survival: float
    n_runs: int
    n_holds: int
    baseline_holds: bool
    cell_level_survival: float | None
    recovered_survival: float | None
    pseudoreplication_flag: bool
    by_choice: list[dict[str, Any]]  # axis, value, n_runs, n_holds, is_baseline
    breakers: list[dict[str, Any]]  # subset of by_choice with failures, sorted
    oat_flips: list[dict[str, Any]]  # one-at-a-time runs that disagree with baseline
    summary: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ClaimResult:
    spec: ClaimSpec
    mode: str  # "run" | "audit"
    verdict: Verdict
    runs: list[RunEvaluation]
    as_labelled: dict[str, Any] | None
    reference_source: str
    data_summary: dict[str, Any]
    grid: GridConfig
    thresholds: Thresholds
    provenance: dict[str, Any]
    methods: str
    total_runtime_s: float
