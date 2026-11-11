"""Provenance: everything needed to reproduce a report, written as JSON."""

from __future__ import annotations

import datetime as dt
import importlib.metadata as md
import json
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from recluster import __version__

PACKAGES = [
    "scanpy",
    "anndata",
    "numpy",
    "scipy",
    "pandas",
    "scikit-learn",
    "igraph",
    "harmonypy",
    "pydeseq2",
    "scikit-image",
    "numba",
    "pydantic",
]


def package_versions() -> dict[str, str]:
    out = {"recluster": __version__}
    for p in PACKAGES:
        try:
            out[p] = md.version(p)
        except md.PackageNotFoundError:
            out[p] = "not installed"
    return out


def environment() -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": package_versions(),
    }


def timestamp() -> str:
    return dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat()


def json_default(o: Any) -> Any:
    """Make numpy / pandas / pydantic objects JSON-serialisable."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, pd.Index | pd.Series):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    if hasattr(o, "model_dump"):
        return o.model_dump()
    if isinstance(o, set):
        return sorted(o)
    raise TypeError(f"not JSON serialisable: {type(o).__name__}")


def dump_json(obj: Any, path: Path) -> None:
    path.write_text(json.dumps(obj, indent=2, default=json_default, allow_nan=False))


def round_floats(obj: Any, ndigits: int = 6) -> Any:
    """Recursively round floats so re-runs compare equal and the JSON stays readable."""
    if isinstance(obj, float):
        if obj != obj:  # NaN
            return None
        return float(f"{obj:.{ndigits}g}")
    if isinstance(obj, dict):
        return {k: round_floats(v, ndigits) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [round_floats(v, ndigits) for v in obj]
    if isinstance(obj, np.floating):
        return round_floats(float(obj), ndigits)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj
