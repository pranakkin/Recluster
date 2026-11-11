"""Doublet detection vs demuxlet ground truth on the Kang 2018 batch-2 lanes.

GEO GSE96583 provides the two 10x lanes of the IFN-beta experiment with *all* barcodes
(including the doublets that the processed object drops) and demuxlet's per-barcode call
(``singlet`` / ``doublet`` / ``ambs``). Recluster scores doublets exactly as the pipeline does
(Scrublet per lane, baseline QC) and the calls are compared.

Caveat stated up front: demuxlet can only see doublets formed by cells from *different*
individuals. With 8 roughly equally mixed donors about 7/8 of real doublets are detectable,
so demuxlet "singlets" contain some same-donor doublets. Agreement is therefore a lower bound.
"""

from __future__ import annotations

import gzip
import hashlib
import urllib.request
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
from anndata import AnnData
from sklearn.metrics import roc_auc_score

from recluster.config import GridConfig
from recluster.pipeline import SAMPLE_COL, Pipeline

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
FILES: dict[str, tuple[str, str]] = {
    "GSE96583_batch2.total.tsne.df.tsv.gz": (
        f"{GEO}/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz",
        "1d57e72e92ca8695250e88cc0f1c3fa8c0be1175d974f8b427c58f1274dc6c09",
    ),
    "GSE96583_batch2.genes.tsv.gz": (
        f"{GEO}/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz",
        "93aa4e9b530ef9d6411ca129b416324c5cc1cc5a01a1fa6ed4f4a845480ed3ca",
    ),
    "GSM2560248_2.1.mtx.gz": (
        f"{GEO}/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_2.1.mtx.gz",
        "32add28a0b3397d9ef3f220b7a6a55e98e60fe7b66fe48d0986d634df8ca0013",
    ),
    "GSM2560248_barcodes.tsv.gz": (
        f"{GEO}/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_barcodes.tsv.gz",
        "d58d8d55cbe4a12757207784b3bc9227bf200c9100ca15131176e9f8159c955e",
    ),
    "GSM2560249_2.2.mtx.gz": (
        f"{GEO}/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_2.2.mtx.gz",
        "8aecc98a7ac4957bbc2570f87ebe8ce97332a5bcdbf557d40bd5aabfd287bdc5",
    ),
    "GSM2560249_barcodes.tsv.gz": (
        f"{GEO}/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_barcodes.tsv.gz",
        "9bb38e080dfae81036fbcd9902c6c6254a4466aa85f70705479be3b2d6679d55",
    ),
}
LANES = {
    "ctrl": ("GSM2560248_2.1.mtx.gz", "GSM2560248_barcodes.tsv.gz"),
    "stim": ("GSM2560249_2.2.mtx.gz", "GSM2560249_barcodes.tsv.gz"),
}


def download_geo(geo_dir: Path) -> None:
    geo_dir.mkdir(parents=True, exist_ok=True)
    for name, (url, sha) in FILES.items():
        path = geo_dir / name
        if not path.exists():
            urllib.request.urlretrieve(url, path)  # noqa: S310 - fixed https URLs
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != sha:
            raise RuntimeError(f"{path}: SHA-256 {digest} != expected {sha}")


def load_batch2(geo_dir: Path) -> AnnData:
    """Both lanes with every barcode, joined to the demuxlet table."""
    download_geo(geo_dir)
    genes = pd.read_csv(geo_dir / "GSE96583_batch2.genes.tsv.gz", sep="\t", header=None)
    calls = pd.read_csv(geo_dir / "GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    calls.index.name = "barcode"
    calls = calls.reset_index()
    parts = []
    for lane, (mtx, bcs) in LANES.items():
        with gzip.open(geo_dir / mtx, "rb") as fh:
            x = scipy.io.mmread(fh).T.tocsr().astype(np.float32)  # cells x genes
        barcodes = pd.read_csv(geo_dir / bcs, sep="\t", header=None)[0].astype(str).to_numpy()
        var = pd.DataFrame(index=pd.Index(genes[1].astype(str).to_numpy()))
        var.index = _unique(var.index)
        a = AnnData(X=x, var=var)
        a.obs_names = [f"{lane}:{b}" for b in barcodes]
        a.obs["lane"] = lane
        a.obs["barcode"] = barcodes
        sub = calls[calls["stim"] == lane].set_index("barcode")
        a.obs["demuxlet"] = sub["multiplets"].reindex(barcodes).to_numpy()
        a.obs["demuxlet_ind"] = sub["ind"].reindex(barcodes).astype("string").to_numpy()
        a.obs["demuxlet_cell"] = sub["cell"].reindex(barcodes).to_numpy()
        parts.append(a)
    out = ad.concat(parts, join="outer", merge="same")
    out.var = out.var[[]]
    return out


def _unique(index: pd.Index) -> pd.Index:
    """Suffix duplicate gene symbols the way anndata does (``NAME-1``, ``NAME-2``)."""
    seen: dict[str, int] = {}
    out = []
    for name in index.astype(str):
        if name in seen:
            seen[name] += 1
            out.append(f"{name}-{seen[name]}")
        else:
            seen[name] = 0
            out.append(name)
    return pd.Index(out)


def doublet_check(geo_dir: Path, *, grid: GridConfig | None = None) -> dict[str, Any]:
    """Scrublet (as run by Recluster, per lane, baseline QC) against demuxlet calls."""
    grid = grid or GridConfig()
    adata = load_batch2(geo_dir)
    pipe = Pipeline(adata, grid, sample=adata.obs["lane"], mito_prefix=None)
    stage = pipe.qc_stage(grid.baseline)
    obs = stage.adata.obs
    out: dict[str, Any] = {
        "n_barcodes": int(adata.n_obs),
        "n_after_qc": int(stage.adata.n_obs),
        "demuxlet_calls_all": adata.obs["demuxlet"].value_counts(dropna=False).to_dict(),
        "warnings": stage.warnings,
        "note": (
            "demuxlet detects only inter-individual doublets (~7/8 of doublets with 8 donors); "
            "'ambs' barcodes are excluded from precision/recall"
        ),
    }
    out["per_lane"] = {}
    for lane in ["all", *sorted(obs[SAMPLE_COL].unique())]:
        sel = obs if lane == "all" else obs[obs[SAMPLE_COL] == lane]
        out["per_lane"][lane] = _agreement(sel)
    return out


def _agreement(obs: pd.DataFrame) -> dict[str, Any]:
    known = obs[obs["demuxlet"].isin(["singlet", "doublet"])]
    truth = (known["demuxlet"] == "doublet").to_numpy()
    pred = known["predicted_doublet"].to_numpy(dtype=bool)
    score = known["doublet_score"].to_numpy(dtype=float)
    tp = int((pred & truth).sum())
    fp = int((pred & ~truth).sum())
    fn = int((~pred & truth).sum())
    tn = int((~pred & ~truth).sum())
    ok = ~np.isnan(score)
    auroc = (
        float(roc_auc_score(truth[ok], score[ok]))
        if truth[ok].any() and (~truth[ok]).any()
        else None
    )
    return {
        "n_cells": int(len(known)),
        "demuxlet_doublet_rate": round(float(truth.mean()), 4),
        "scrublet_called_rate": round(float(pred.mean()), 4),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": round(tp / (tp + fp), 4) if tp + fp else None,
        "recall": round(tp / (tp + fn), 4) if tp + fn else None,
        "auroc_of_score": None if auroc is None else round(auroc, 4),
        "ambiguous_excluded": int((obs["demuxlet"] == "ambs").sum()),
    }


def summarize_for_markdown(d: dict[str, Any]) -> str:
    rows = [
        "| Lane | Cells | demuxlet doublet rate | Scrublet called | Precision | Recall | AUROC |",
        "|---|---|---|---|---|---|---|",
    ]
    for lane, a in d["per_lane"].items():
        rows.append(
            f"| {lane} | {a['n_cells']:,} | {a['demuxlet_doublet_rate']:.1%} | "
            f"{a['scrublet_called_rate']:.1%} | "
            f"{_f(a['precision'])} | {_f(a['recall'])} | {_f(a['auroc_of_score'])} |"
        )
    return "\n".join(rows) + "\n"


def _f(x: float | None) -> str:
    return "–" if x is None else f"{x:.2f}"
