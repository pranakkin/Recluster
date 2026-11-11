"""One deterministic pipeline run, from raw counts to cluster labels.

QC -> (doublet scoring, cached per QC setting) -> normalisation -> HVG -> PCA ->
(harmony) -> kNN -> leiden. Only cluster assignments are kept per run; expression used by the
evaluators is computed once, elsewhere, from the raw counts.
"""

from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
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
    effective: dict[str, object] = field(default_factory=dict)


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
    doublet_scores: pd.Series  # Scrublet score per kept cell (NaN if scoring failed)
    warnings: list[str] = field(default_factory=list)
    effective: dict[str, object] = field(default_factory=dict)

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
            # float64, so a rounded threshold prints as 5.3 rather than 5.300000190734863
            mt = np.asarray(x[:, is_mt].sum(axis=1)).ravel().astype(np.float64)
            total64 = total.astype(np.float64)
            with np.errstate(divide="ignore", invalid="ignore"):
                a.obs["pct_counts_mt"] = np.where(total64 > 0, 100.0 * mt / total64, 0.0)
        else:
            a.obs["pct_counts_mt"] = 0.0

    def effective_qc(self, params: PipelineParams) -> tuple[int, float, dict[str, object]]:
        """QC thresholds after the cell-loss cap.

        A threshold that alone would remove more than ``qc_max_cell_loss`` of the cells is not
        a defensible alternative on this dataset; it is pulled back to the value that removes
        exactly that fraction, and the run records both values.
        """
        cap = self.cfg.qc_max_cell_loss
        obs = self.counts.obs
        n_genes = obs["n_genes_by_counts"].to_numpy()
        pct_mt = obs["pct_counts_mt"].to_numpy()
        capped: dict[str, object] = {}
        min_genes = int(params.min_genes)
        if n_genes.size and (n_genes < min_genes).mean() > cap:
            min_genes = int(np.quantile(n_genes, cap))
            capped["min_genes"] = min_genes
        max_mito = float(params.max_mito_pct)
        if self.mito_prefix is not None and pct_mt.size and (pct_mt > max_mito).mean() > cap:
            max_mito = round(float(np.quantile(pct_mt.astype(np.float64), 1 - cap)), 1)
            capped["max_mito_pct"] = max_mito
        return min_genes, max_mito, capped

    def effective_params(self, params: PipelineParams) -> PipelineParams:
        """``params`` with ``min_genes`` / ``max_mito_pct`` replaced by their capped values.

        Grid points whose thresholds are pulled back to the same values run the same pipeline;
        the engine uses this to drop such duplicates.
        """
        _, _, capped = self.effective_qc(params)
        return params.model_copy(update=capped) if capped else params

    def qc_stage(self, params: PipelineParams) -> QCStage:
        """The QC stage for ``params``, with warnings and cap info composed for these params.

        Stages are cached by effective thresholds, so grid points that cap to the same values
        share one. The cached object holds only what is true for all of them; the lines naming
        this call's requested values, and which of its axes were capped (a request can equal
        the capped value exactly), go on a shallow copy returned here.
        """
        stage, capped = self._qc_cached(params)
        return replace(
            stage,
            warnings=stage.warnings + self._cap_warnings(params, capped),
            effective={**stage.effective, "capped": capped},
        )

    def _qc_cached(self, params: PipelineParams) -> tuple[QCStage, dict[str, object]]:
        """Cached stage for the effective thresholds of ``params``, and which axes were capped."""
        min_genes, max_mito, capped = self.effective_qc(params)
        mito_key = max_mito if self.mito_prefix is not None else float("inf")
        key = (min_genes, mito_key)
        stage = self._qc_cache.get(key)
        if stage is None:
            stage = self._build_qc_stage(min_genes, max_mito, capped)
            self._qc_cache[key] = stage
        return stage, capped

    def _build_qc_stage(
        self, min_genes: int, max_mito: float, capped: dict[str, object]
    ) -> QCStage:
        """Filter cells and genes at the given thresholds and score doublets."""
        t0 = time.perf_counter()
        warns: list[str] = []
        obs = self.counts.obs
        keep = (obs["n_genes_by_counts"] >= min_genes) & (obs["pct_counts_mt"] <= max_mito)
        a = self.counts[keep.values].copy()
        sc.pp.filter_genes(a, min_cells=3)
        if self.mito_prefix is None:
            warns.append("no mitochondrial genes found; max_mito_pct has no effect")
        fallback_info = self._score_doublets(a, warns)
        return QCStage(
            adata=a,
            n_input=int(self.counts.n_obs),
            n_removed_qc=int(self.counts.n_obs - a.n_obs),
            runtime_s=time.perf_counter() - t0,
            warnings=warns,
            effective={
                "min_genes": min_genes,
                "max_mito_pct": max_mito,
                "capped": capped,
                "scrublet_fallback_subsets": len(fallback_info),
            },
        )

    def _cap_warnings(self, params: PipelineParams, capped: dict[str, object]) -> list[str]:
        """One line per capped axis, naming the value ``params`` asked for."""
        return [
            f"{axis}={getattr(params, axis)} would remove more than "
            f"{self.cfg.qc_max_cell_loss:.0%} of cells; capped to {value}"
            for axis, value in capped.items()
        ]

    def _score_doublets(self, a: AnnData, warns: list[str]) -> list[dict[str, object]]:
        """Scrublet per sample when every sample has >= 100 cells; otherwise all cells together.

        Each sample is scored in its own call (doublets form within a 10x lane, not across
        lanes). scanpy's ``batch_key`` loop is avoided on purpose: with scipy >= 1.17 it fails
        on sparse arrays. Returns one record per subset where the threshold fallback was used.
        """
        groups: np.ndarray | None = None
        if self.has_sample:
            sizes = a.obs[SAMPLE_COL].value_counts()
            if (sizes < 100).any():
                warns.append(
                    f"{int((sizes < 100).sum())} sample(s) have <100 cells; "
                    "scrublet run on all cells together"
                )
            else:
                groups = a.obs[SAMPLE_COL].astype(str).to_numpy()
        scores = np.full(a.n_obs, np.nan)
        pred = np.zeros(a.n_obs, dtype=bool)
        subsets = (
            [np.arange(a.n_obs)]
            if groups is None
            else [np.flatnonzero(groups == g) for g in np.unique(groups)]
        )
        n_failed = 0
        fallbacks: list[dict[str, object]] = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for idx in subsets:
                try:
                    s, p, info = self._scrublet_subset(a, idx)
                except Exception as exc:  # noqa: BLE001 - keep going, report in provenance
                    n_failed += 1
                    warns.append(f"scrublet failed on one subset ({exc.__class__.__name__}: {exc})")
                    continue
                scores[idx], pred[idx] = s, p
                if info.get("fallback"):
                    fallbacks.append(info)
        if fallbacks:
            warns.append(
                f"scrublet's automatic threshold called almost no doublets in {len(fallbacks)} "
                f"subset(s); the top expected fraction by score was called instead "
                f"(expected rate {100 * self.cfg.scrublet_expected_rate_per_1000:g}% "
                "per 1,000 cells)"
            )
        a.obs["doublet_score"] = scores
        a.obs["predicted_doublet"] = pred
        return fallbacks

    def _scrublet_subset(
        self, a: AnnData, idx: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
        """Scrublet on one subset, with a rate-based fallback for its automatic threshold.

        Scrublet's threshold comes from the valley of the simulated-doublet score histogram and
        can land above every observed score, calling nothing. When the called fraction is below
        ``scrublet_fallback_trigger`` times the expected doublet rate (0.8 % per 1,000 cells
        recovered, the usual 10x rule of thumb), the top expected fraction by score is called.
        The expectation assumes each scored subset was captured on its own lane; for samples
        demultiplexed from a pooled lane the true rate is the lane's.
        """
        sub = AnnData(
            X=sp.csr_matrix(a.X[idx]), obs=a.obs.iloc[idx][[]].copy(), var=a.var[[]].copy()
        )
        sc.pp.scrublet(
            sub,
            random_state=self.cfg.scrublet_seed,
            use_approx_neighbors=False,  # exact kNN: deterministic, and avoids a scanpy/scipy bug
            verbose=False,
        )
        scores = sub.obs["doublet_score"].to_numpy(dtype=float)
        pred_raw = sub.obs["predicted_doublet"]
        pred = pred_raw.fillna(False).astype(bool).to_numpy()
        n = len(idx)
        expected = min(0.30, self.cfg.scrublet_expected_rate_per_1000 * n / 1000.0)
        called = float(pred.mean()) if n else 0.0
        info: dict[str, object] = {
            "n_cells": n,
            "expected_rate": expected,
            "called_fraction_auto": called,
            "fallback": False,
        }
        if (
            n
            and expected > 0
            and (pred_raw.isna().any() or called < self.cfg.scrublet_fallback_trigger * expected)
        ):
            k = int(round(expected * n))
            order = np.argsort(-np.nan_to_num(scores, nan=-1.0), kind="stable")
            pred = np.zeros(n, dtype=bool)
            pred[order[:k]] = True
            info.update(
                {"fallback": True, "k": k, "kth_score": float(scores[order[k - 1]]) if k else None}
            )
        return scores, pred, info

    # ----------------------------------------------------------------- run

    def run(self, params: PipelineParams) -> RunResult:
        t0 = time.perf_counter()
        stage, capped = self._qc_cached(params)
        warns = stage.warnings + self._cap_warnings(params, capped)
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
        effective = {**stage.effective, "capped": capped}
        return RunResult(
            params=params,
            clusters=clusters,
            n_input=stage.n_input,
            n_removed_qc=stage.n_removed_qc,
            n_removed_doublet=n_removed_doublet,
            n_clusters=int(clusters.nunique()),
            runtime_s=time.perf_counter() - t0,
            doublet_scores=a.obs["doublet_score"].astype(float),
            warnings=warns,
            effective=effective,
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
