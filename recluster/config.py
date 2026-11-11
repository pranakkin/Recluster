"""Grid configuration and verdict thresholds. Everything here is explicit and overridable."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field

Normalization = Literal["log_total", "pearson"]


class PipelineParams(BaseModel):
    """One point in the perturbation grid: every choice a single pipeline run makes."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    min_genes: int = 200
    max_mito_pct: float = 10.0
    doublet_removal: bool = True
    normalization: Normalization = "log_total"
    n_hvgs: int = 2000
    n_pcs: int = 30
    resolution: float = 1.0
    batch_correction: bool = True  # forced to False by build_grid when there is no batch column
    seed: int = 0

    def diff(self, other: PipelineParams) -> dict[str, object]:
        """Fields where ``self`` differs from ``other`` (used for 'what broke it')."""
        return {k: v for k, v in self.model_dump().items() if v != getattr(other, k)}

    def label(self, baseline: PipelineParams) -> str:
        d = self.diff(baseline)
        if not d:
            return "baseline"
        return ", ".join(f"{k}={v}" for k, v in d.items())


class GridAxes(BaseModel):
    """Defensible alternatives for each pipeline choice. Baseline values need not be listed."""

    model_config = ConfigDict(extra="forbid")

    min_genes: list[int] = [100, 200, 500]
    max_mito_pct: list[float] = [5.0, 10.0, 20.0]
    doublet_removal: list[bool] = [False, True]
    normalization: list[Normalization] = ["log_total", "pearson"]
    n_hvgs: list[int] = [1000, 2000, 3000]
    n_pcs: list[int] = [20, 30, 50]
    resolution: list[float] = [0.3, 0.5, 1.0, 1.5, 2.0]
    batch_correction: list[bool] = [False, True]
    seed: list[int] = [0, 1, 2]


class GridConfig(BaseModel):
    """The bounded experimental design plus fixed pipeline settings."""

    model_config = ConfigDict(extra="forbid")

    baseline: PipelineParams = Field(default_factory=PipelineParams)
    axes: GridAxes = Field(default_factory=GridAxes)
    n_random: int = Field(
        default=6, ge=0, description="random full-factorial draws added to one-at-a-time"
    )
    random_seed: int = 0
    n_neighbors: int = 15
    mito_prefix: str = Field(default="auto", description="'auto' tries MT- then mt-")
    scrublet_seed: int = Field(default=0, description="fixed so the doublet axis is purely on/off")
    max_runs: int | None = Field(default=None, description="hard cap on grid size, for quick runs")
    n_cpus: int = Field(default=4, ge=1, description="threads for pydeseq2")

    @classmethod
    def from_yaml(cls, path: str | Path) -> GridConfig:
        with open(path) as fh:
            return cls.model_validate(yaml.safe_load(fh) or {})


class Thresholds(BaseModel):
    """Every cutoff used to turn computed numbers into a verdict. Printed in each report."""

    model_config = ConfigDict(extra="forbid")

    robust_min_survival: float = Field(default=0.80, ge=0, le=1)
    fragile_min_survival: float = Field(default=0.20, ge=0, le=1)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    marker_auroc_min: float = Field(default=0.80, ge=0.5, le=1)
    jaccard_min: float = Field(default=0.50, ge=0, le=1)
    label_precision_min: float = Field(default=0.50, ge=0, le=1)
    min_samples_present: int = Field(default=2, ge=1)
    min_cells_per_sample: int = Field(default=5, ge=1)
    min_cluster_size: int = Field(
        default=10, ge=1, description="clusters smaller than this are never candidates"
    )
    composition_test: Literal["logit_t", "mannwhitney"] = "logit_t"
    pseudobulk_min_total_counts: int = Field(default=10, ge=0)
    pseudoreplication_cell_survival: float = Field(default=0.80, ge=0, le=1)
    pseudoreplication_replicate_survival: float = Field(default=0.50, ge=0, le=1)

    @classmethod
    def from_yaml(cls, path: str | Path) -> Thresholds:
        with open(path) as fh:
            return cls.model_validate(yaml.safe_load(fh) or {})
