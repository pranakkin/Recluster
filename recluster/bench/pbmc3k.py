"""PBMC3k loading and a marker-named reference clustering used to build planted errors.

The reference labels are produced by Recluster's own baseline pipeline at resolution 1.0 and
named by the standard PBMC marker sets (as in the scanpy tutorial). They are only used to
*construct* planted errors; no benchmark number depends on them being perfect.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from recluster.config import GridConfig, PipelineParams
from recluster.data import detect_mito_prefix, extract_counts
from recluster.matching import lognorm_matrix, marker_score
from recluster.pipeline import Pipeline

MARKERS: dict[str, list[str]] = {
    "CD4 T": ["IL7R", "CD3D", "CCR7"],
    "CD8 T": ["CD8A", "CD8B", "CD3D"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7"],
    "DC": ["FCER1A", "CST3", "CLEC10A"],
    "Mk": ["PPBP", "PF4"],
}


def load_pbmc3k(data_dir: str | Path = "data") -> AnnData:
    """Raw PBMC3k counts (downloaded by scanpy into ``data_dir`` on first use)."""
    sc.settings.datasetdir = str(data_dir)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    return adata


def reference_labels(adata: AnnData, *, resolution: float = 1.0, seed: int = 0) -> pd.Series:
    """Cluster with the baseline pipeline and name clusters by the strongest marker set.

    Cells removed by QC/doublet filtering get the label ``"filtered"``.
    """
    counts, _ = extract_counts(adata)
    pipe = Pipeline(counts, GridConfig(), mito_prefix=detect_mito_prefix(counts.var_names))
    run = pipe.run(PipelineParams(resolution=resolution, seed=seed, batch_correction=False))
    genes = sorted({g for gs in MARKERS.values() for g in gs if g in counts.var_names})
    ln = lognorm_matrix(counts, genes)
    scores = pd.DataFrame(
        {t: marker_score(ln, [g for g in gs if g in genes], []) for t, gs in MARKERS.items()}
    )
    per_cluster = scores.loc[run.clusters.index].groupby(run.clusters.to_numpy()).mean()
    names = per_cluster.idxmax(axis=1)
    labels = run.clusters.map(names).astype(str)
    out = pd.Series("filtered", index=adata.obs_names, name="reference")
    out.loc[labels.index] = labels.to_numpy()
    return out


def as_processed(
    adata: AnnData, labels: pd.Series | None = None, label_col: str = "tool_label"
) -> AnnData:
    """Return an object shaped like another tool's output: log-normalised .X, counts in layers."""
    out = adata.copy()
    out.layers["counts"] = out.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(out, target_sum=1e4)
        sc.pp.log1p(out)
    if labels is not None:
        out.obs[label_col] = labels.loc[out.obs_names].astype(str).to_numpy()
    return out


def expressed_genes(
    adata: AnnData, min_fraction: float, rng: np.random.Generator, n: int
) -> np.ndarray:
    """Sample ``n`` genes expressed in at least ``min_fraction`` of cells."""
    frac = np.asarray((adata.X > 0).mean(axis=0)).ravel()
    pool = np.flatnonzero(frac >= min_fraction)
    return rng.choice(pool, size=min(n, len(pool)), replace=False)
