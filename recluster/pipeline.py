"""One deterministic pipeline run, from raw counts to cluster labels.

QC -> (doublet scoring, cached per QC setting) -> normalisation -> HVG -> PCA ->
(harmony) -> kNN -> leiden. Only cluster assignments are kept per run; expression used by the
evaluators is computed once, elsewhere, from the raw counts.
"""

from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from recluster.config import GridConfig, PipelineParams

log = logging.getLogger(__name__)

SAMPLE_COL = "_recluster_sample"
BATCH_COL = "_recluster_batch"
CLUSTER_COL = "cluster"


@dataclass
class QCStage:
    """Counts after cell/gene QC with doublet scores attached."""

    adata: AnnData
    n_input: int
    n_removed_qc: int
    runtime_s: float
    warnings: list[str] = field(default_factory=list)


@dataclass
class RunResult:
    """Everything a grid run contributes downstream."""

    params: PipelineParams
    clusters: pd.Series  # index: kept barcodes; values: cluster id as str
    n_input: int
    n_removed_qc: int
    n_removed_doublet: int
    n_clusters: int
    runtime_s: float
    warnings: list[str] = field(default_factory=list)

    @property
    def n_cells(self) -> int:
        return int(len(self.clusters))


class Pipeline:
    """Holds the raw counts and a QC-stage cache; runs grid points."""

    def __init__(
        self,
        counts: AnnData,
        cfg: GridConfig,
        *,
        sample: pd.Series | None = None,
        batch: pd.Series | None = None,
        mito_prefix: str | None = None,
    ) -> None:
        self.counts = counts
        self.cfg = cfg
        self.mito_prefix = mito_prefix
        if sample is not None:
            counts.obs[SAMPLE_COL] = sample.loc[counts.obs_names].astype(str).values
        if batch is not None:
            counts.obs[BATCH_COL] = batch.loc[counts.obs_names].astype(str).values
        self.has_sample = sample is not None
        self.has_batch = batch is not None
        self._qc_cache: dict[tuple[int, float], QCStage] = {}
        self._precompute_qc_metrics()

    # ------------------------------------------------------------------ QC

    def _precompute_qc_metrics(self) -> None:
        a = self.counts
        x = a.X
        a.obs["n_genes_by_counts"] = np.asarray((x > 0).sum(axis=1)).ravel()
        total = np.asarray(x.sum(axis=1)).ravel()
        a.obs["total_counts"] = total
        if self.mito_prefix is not None:
            is_mt = np.char.startswith(np.asarray(a.var_names, dtype=str), self.mito_prefix)
            mt = np.asarray(x[:, is_mt].sum(axis=1)).ravel()
            with np.errstate(divide="ignore", invalid="ignore"):
                a.obs["pct_counts_mt"] = np.where(total > 0, 100.0 * mt / total, 0.0)
        else:
            a.obs["pct_counts_mt"] = 0.0

    def qc_stage(self, params: PipelineParams) -> QCStage:
        key = (int(params.min_genes), float(params.max_mito_pct))
        if key in self._qc_cache:
            return self._qc_cache[key]
        t0 = time.perf_counter()
        warns: list[str] = []
        obs = self.counts.obs
        keep = (obs["n_genes_by_counts"] >= params.min_genes) & (
            obs["pct_counts_mt"] <= params.max_mito_pct
        )
        a = self.counts[keep.values].copy()
        sc.pp.filter_genes(a, min_cells=3)
        if self.mito_prefix is None:
            warns.append("no mitochondrial genes found; max_mito_pct has no effect")
        self._score_doublets(a, warns)
        stage = QCStage(
            adata=a,
            n_input=int(self.counts.n_obs),
            n_removed_qc=int(self.counts.n_obs - a.n_obs),
            runtime_s=time.perf_counter() - t0,
            warnings=warns,
        )
        self._qc_cache[key] = stage
        return stage

    def _score_doublets(self, a: AnnData, warns: list[str]) -> None:
        """Scrublet per sample when possible; falls back to all cells together."""
        batch_key = SAMPLE_COL if self.has_sample else None
        if batch_key is not None:
            sizes = a.obs[batch_key].value_counts()
            if (sizes < 100).any():
                warns.append(
                    f"{int((sizes < 100).sum())} sample(s) have <100 cells; "
                    "scrublet run on all cells together"
                )
                batch_key = None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                sc.pp.scrublet(
                    a, batch_key=batch_key, random_state=self.cfg.scrublet_seed, verbose=False
                )
            except Exception as exc:  # noqa: BLE001 - we want to keep going and report
                warns.append(
                    f"scrublet failed ({exc.__class__.__name__}: {exc}); no doublets called"
                )
                a.obs["doublet_score"] = np.nan
                a.obs["predicted_doublet"] = False
        pred = a.obs["predicted_doublet"]
        if pred.isna().any():
            warns.append("scrublet could not set a threshold for some cells; treated as singlets")
        a.obs["predicted_doublet"] = pred.fillna(False).astype(bool)

    # ----------------------------------------------------------------- run

    def run(self, params: PipelineParams) -> RunResult:
        t0 = time.perf_counter()
        stage = self.qc_stage(params)
        warns = list(stage.warnings)
        a = stage.adata
        if params.doublet_removal:
            keep = ~a.obs["predicted_doublet"].values
            n_removed_doublet = int((~keep).sum())
            a = a[keep].copy()
        else:
            n_removed_doublet = 0
            a = a.copy()

        n_hvgs = min(params.n_hvgs, a.n_vars)
        if n_hvgs < params.n_hvgs:
            warns.append(f"only {a.n_vars} genes after QC; n_hvgs reduced to {n_hvgs}")
        n_pcs = min(params.n_pcs, n_hvgs - 1, a.n_obs - 1)
        if n_pcs < params.n_pcs:
            warns.append(f"n_pcs reduced to {n_pcs} (limited by genes/cells)")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if params.normalization == "log_total":
                sc.pp.normalize_total(a, target_sum=1e4)
                sc.pp.log1p(a)
                sc.pp.highly_variable_genes(a, n_top_genes=n_hvgs, flavor="seurat")
                a = a[:, a.var["highly_variable"].values].copy()
                sc.pp.scale(a, max_value=10)
            elif params.normalization == "pearson":
                sc.experimental.pp.highly_variable_genes(
                    a, flavor="pearson_residuals", n_top_genes=n_hvgs
                )
                a = a[:, a.var["highly_variable"].values].copy()
                sc.experimental.pp.normalize_pearson_residuals(a)
            else:  # pragma: no cover - pydantic restricts the literal
                raise ValueError(params.normalization)
            sc.pp.pca(a, n_comps=n_pcs, svd_solver="arpack", random_state=params.seed)

            rep = "X_pca"
            if params.batch_correction and self.has_batch:
                rep = self._harmony(a, params.seed, warns)

            sc.pp.neighbors(
                a,
                n_neighbors=self.cfg.n_neighbors,
                use_rep=rep,
                random_state=params.seed,
                transformer="sklearn",
            )
            sc.tl.leiden(
                a,
                resolution=params.resolution,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=params.seed,
                key_added=CLUSTER_COL,
            )
        clusters = a.obs[CLUSTER_COL].astype(str)
        clusters.name = CLUSTER_COL
        return RunResult(
            params=params,
            clusters=clusters,
            n_input=stage.n_input,
            n_removed_qc=stage.n_removed_qc,
            n_removed_doublet=n_removed_doublet,
            n_clusters=int(clusters.nunique()),
            runtime_s=time.perf_counter() - t0,
            warnings=warns,
        )

    @staticmethod
    def _harmony(a: AnnData, seed: int, warns: list[str]) -> str:
        import harmonypy  # local import: optional heavy dependency path

        if a.obs[BATCH_COL].nunique() < 2:
            warns.append("batch column has a single level; harmony skipped")
            return "X_pca"
        ho = harmonypy.run_harmony(
            a.obsm["X_pca"], a.obs[[BATCH_COL]], BATCH_COL, random_state=seed, verbose=False
        )
        z = ho.Z_corr
        z = z.cpu().numpy() if hasattr(z, "cpu") else np.asarray(z)
        if z.shape[0] != a.n_obs:
            z = z.T
        a.obsm["X_harmony"] = np.ascontiguousarray(z, dtype=np.float32)
        return "X_harmony"
