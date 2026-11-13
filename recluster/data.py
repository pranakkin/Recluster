"""Loading AnnData, locating raw counts, hashing inputs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from anndata import AnnData


class DataError(ValueError):
    """Raised when the input AnnData cannot be used."""


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """Streaming SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


AUDIT_COUNTS_IN_X = (
    "audit mode expects a processed object whose raw counts live in layers['counts'] "
    "or .raw; "
    ".X itself looks like raw counts, so use `recluster run` instead"
)


def load_adata(path: str | Path) -> AnnData:
    path = Path(path)
    if not path.exists():
        raise DataError(f"{path} does not exist")
    try:
        return ad.read_h5ad(path)
    except OSError as exc:
        raise DataError(
            f"{path} could not be read as an AnnData .h5ad file ({exc}). "
            "--data takes the .h5ad and --claim the YAML."
        ) from exc


def _looks_like_counts(x: sp.spmatrix | np.ndarray, n_check: int = 2000) -> bool:
    """True if a sample of the non-zero values are non-negative integers."""
    if sp.issparse(x):
        data = np.asarray(x.data)
    else:
        data = np.asarray(x).ravel()
        data = data[data != 0]
    if data.size == 0:
        return False
    sample = data[: min(n_check, data.size)]
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def extract_counts(adata: AnnData) -> tuple[AnnData, str]:
    """Return a new AnnData whose ``.X`` is raw counts (CSR) for all cells.

    Search order: ``layers['counts']``, ``.raw.X``, ``.X``. Returns the source used.
    The returned object carries ``.obs`` (copied) and gene names of the counts matrix.
    """
    if "counts" in adata.layers and _looks_like_counts(adata.layers["counts"]):
        x, var = adata.layers["counts"], adata.var.copy()
        source = "layers['counts']"
    elif adata.raw is not None and _looks_like_counts(adata.raw.X):
        x, var = adata.raw.X, adata.raw.var.copy()
        source = ".raw.X"
    elif _looks_like_counts(adata.X):
        x, var = adata.X, adata.var.copy()
        source = ".X"
    else:
        raise DataError(
            "no raw count matrix found: checked layers['counts'], .raw.X and .X; "
            "all contain non-integer values. Recluster needs raw counts to re-derive results."
        )
    x = sp.csr_matrix(x) if not sp.isspmatrix_csr(x) else x.copy()
    x = x.astype(np.float32)
    counts = AnnData(X=x, obs=adata.obs.copy(), var=var[[]].copy())
    counts.var_names_make_unique()
    return counts, source


def detect_mito_prefix(var_names: list[str] | np.ndarray, requested: str = "auto") -> str | None:
    """Pick the mitochondrial gene prefix; None if no mito genes are found."""
    candidates = ["MT-", "mt-"] if requested == "auto" else [requested]
    names = np.asarray(var_names, dtype=str)
    for c in candidates:
        if np.char.startswith(names, c).any():
            return c
    return None
