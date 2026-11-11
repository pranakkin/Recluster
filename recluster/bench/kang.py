"""Kang et al. 2018 (Nat Biotechnol 36:89): PBMCs from 8 lupus patients, IFN-beta vs control.

Source: the processed object distributed with pertpy (``kang_2018.h5ad`` on
exampledata.scverse.org),
raw integer counts in .X, columns ``label`` (ctrl/stim), ``replicate`` (8 donors), ``cell_type``.
This object contains singlets only; the demuxlet multiplet calls from GEO are *not* in it.
"""

from __future__ import annotations

import hashlib
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from recluster.spec import ClaimSpec

KANG_URL = "https://exampledata.scverse.org/pertpy/kang_2018.h5ad"
KANG_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"
COLS = {"condition": "label", "sample": ["replicate", "label"], "labels": "cell_type"}
REQUIRED = {"label": {"ctrl", "stim"}, "replicate": None, "cell_type": None}

MONO = {"markers": ["CD14", "LYZ", "S100A8"], "label_value": "CD14+ Monocytes"}
CD4T = {"markers": ["IL7R", "CD3D", "CCR7"], "label_value": "CD4 T cells"}
NK = {"markers": ["GNLY", "NKG7", "KLRD1"], "label_value": "NK cells"}
ISGS = ["ISG15", "IFI6", "IFIT3", "MX1"]


@dataclass
class KangCase:
    key: str
    title: str
    how_built: str
    data_path: Path
    spec: ClaimSpec
    expected: str
    passes_rule: str
    passes: Any  # callable(result) -> bool
    extra: dict[str, Any] = field(default_factory=dict)
    role: str = "true"  # true | null | false


def download_kang(path: Path) -> Path:
    """Download the h5ad if missing and verify its SHA-256."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(KANG_URL, path)  # noqa: S310 - fixed https URL
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != KANG_SHA256:
        raise RuntimeError(f"{path}: SHA-256 {digest} does not match the expected {KANG_SHA256}")
    return path


def check_kang(adata: AnnData) -> dict[str, Any]:
    """Verify the columns the benchmark relies on and report what is (not) there."""
    problems = []
    for col, values in REQUIRED.items():
        if col not in adata.obs.columns:
            problems.append(f"missing obs column {col!r}")
        elif values is not None and set(adata.obs[col].astype(str)) != values:
            problems.append(
                f"{col!r} has values {sorted(set(adata.obs[col].astype(str)))}, "
                f"expected {sorted(values)}"
            )
    if problems:
        raise RuntimeError(
            "Kang 2018 object does not have the expected columns: " + "; ".join(problems)
        )
    doublet_cols = [
        c for c in adata.obs.columns if "multiplet" in c.lower() or "doublet" in c.lower()
    ]
    mito = [g for g in adata.var_names if str(g).upper().startswith("MT-")]
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "donors": sorted(adata.obs["replicate"].astype(str).unique().tolist()),
        "conditions": sorted(adata.obs["label"].astype(str).unique().tolist()),
        "cell_types": adata.obs["cell_type"].astype(str).value_counts().to_dict(),
        "doublet_columns": doublet_cols,
        "n_mito_genes": len(mito),
    }


# ------------------------------------------------------------- constructions


def fake_split_within_control(adata: AnnData, *, seed: int = 0) -> AnnData:
    """Control cells only; donors split into two fake groups A/B (4 vs 4) by seeded shuffle."""
    rng = np.random.default_rng(seed)
    ctrl = adata[(adata.obs["label"].astype(str) == "ctrl").to_numpy()].copy()
    donors = sorted(ctrl.obs["replicate"].astype(str).unique())
    rng.shuffle(donors)
    group = {d: ("A" if i < len(donors) // 2 else "B") for i, d in enumerate(donors)}
    ctrl.obs["fake_group"] = ctrl.obs["replicate"].astype(str).map(group).astype("category")
    ctrl.uns["fake_split"] = {
        "A": sorted(d for d, g in group.items() if g == "A"),
        "B": sorted(d for d, g in group.items() if g == "B"),
    }
    return ctrl


def donor_level_shuffle(adata: AnnData, *, seed: int = 0) -> AnnData:
    """All cells; each of the 16 donor x condition samples gets a fake label A/B, balanced so
    each fake group holds 4 ctrl and 4 stim samples."""
    rng = np.random.default_rng(seed)
    out = adata.copy()
    donors = sorted(out.obs["replicate"].astype(str).unique())
    assign: dict[tuple[str, str], str] = {}
    for cond in ("ctrl", "stim"):
        order = list(donors)
        rng.shuffle(order)
        for i, d in enumerate(order):
            assign[(d, cond)] = "A" if i < len(order) // 2 else "B"
    keys = list(zip(out.obs["replicate"].astype(str), out.obs["label"].astype(str), strict=True))
    out.obs["fake_group"] = pd.Categorical([assign[k] for k in keys])
    out.uns["fake_shuffle"] = {f"{d}|{c}": g for (d, c), g in assign.items()}
    return out


# ------------------------------------------------------------------- claims


def _spec(**kw: Any) -> ClaimSpec:
    return ClaimSpec.model_validate(kw)


def _de(
    name: str,
    target: dict[str, Any],
    genes: list[str],
    direction: str,
    cols: dict[str, Any],
    a: str,
    b: str,
    require: str = "all",
) -> ClaimSpec:
    return _spec(
        claim_type="differential_expression",
        name=name,
        columns=cols,
        target=target,
        comparison={"group_a": a, "group_b": b, "direction": direction},
        genes=genes,
        require=require,
    )


def kang_cases(
    real_path: Path, split_path: Path, shuffle_path: Path, split_groups: dict[str, list[str]]
) -> list[KangCase]:
    """All Kang claims with their expectations fixed before running."""
    fake_cols = {
        "condition": "fake_group",
        "sample": ["replicate", "fake_group"],
        "labels": "cell_type",
    }
    shuffle_cols = {
        "condition": "fake_group",
        "sample": ["replicate", "label"],
        "labels": "cell_type",
    }
    robust = ("verdict is ROBUST", lambda r: r.verdict.label == "ROBUST")
    not_supported = ("verdict is NOT SUPPORTED", lambda r: r.verdict.label == "NOT SUPPORTED")
    real = (
        "Unmodified Kang 2018 object; condition = label (stim vs ctrl), "
        "replicate = donor x condition (16 samples)."
    )
    split = (
        "Control cells only. Donors shuffled with seed 0 and split into fake groups "
        f"A = {', '.join(split_groups['A'])} and B = {', '.join(split_groups['B'])}; "
        "replicate = donor. "
        "No treatment differs between A and B, so any claimed difference is false by construction."
    )
    shuffle = (
        "All cells. Each of the 16 donor x condition samples was assigned a fake group so that "
        "A and B each contain 4 ctrl and 4 stim samples (seed 0). Condition is balanced across "
        "the fake groups, "
        "so claimed differences between A and B are false by construction."
    )
    cases = [
        KangCase(
            "kang_isg15_mono",
            "ISG15 up in CD14+ monocytes after IFN-beta",
            real,
            real_path,
            _de(
                "ISG15 up in CD14+ monocytes, stim vs ctrl (Kang 2018)",
                MONO,
                ["ISG15"],
                "up_in_a",
                COLS,
                "stim",
                "ctrl",
            ),
            "ROBUST",
            *robust,
        ),
        KangCase(
            "kang_isgs_mono",
            "Four interferon-stimulated genes up in CD14+ monocytes",
            real,
            real_path,
            _de(
                "ISG15, IFI6, IFIT3 and MX1 all up in CD14+ monocytes, stim vs ctrl (Kang 2018)",
                MONO,
                ISGS,
                "up_in_a",
                COLS,
                "stim",
                "ctrl",
            ),
            "ROBUST",
            *robust,
        ),
        KangCase(
            "kang_isg15_cd4t",
            "ISG15 up in CD4 T cells after IFN-beta",
            real,
            real_path,
            _de(
                "ISG15 up in CD4 T cells, stim vs ctrl (Kang 2018)",
                CD4T,
                ["ISG15"],
                "up_in_a",
                COLS,
                "stim",
                "ctrl",
            ),
            "ROBUST",
            *robust,
        ),
        KangCase(
            "kang_mono_exists",
            "CD14+ monocytes exist across donors",
            real,
            real_path,
            _spec(
                claim_type="population_exists",
                name="CD14+ monocytes exist (Kang 2018)",
                columns={"sample": ["replicate", "label"], "labels": "cell_type"},
                target={**MONO, "min_cells": 20},
            ),
            "ROBUST",
            *robust,
        ),
        KangCase(
            "kang_nk_exists",
            "NK cells exist across donors",
            real,
            real_path,
            _spec(
                claim_type="population_exists",
                name="NK cells exist (Kang 2018)",
                columns={"sample": ["replicate", "label"], "labels": "cell_type"},
                target={**NK, "negative_markers": ["CD3D"], "min_cells": 20},
            ),
            "ROBUST",
            *robust,
        ),
        KangCase(
            "kang_mono_composition",
            "CD14+ monocyte proportion differs stim vs ctrl",
            real
            + " No compositional change is planted or known: the same donors' PBMCs were split "
            "and one half stimulated for 6 h.",
            real_path,
            _spec(
                claim_type="composition_shift",
                name="CD14+ monocytes differ in abundance, stim vs ctrl (Kang 2018)",
                columns=COLS,
                target=MONO,
                comparison={"group_a": "stim", "group_b": "ctrl", "direction": "different"},
            ),
            "NOT SUPPORTED (no known shift; reported as measured)",
            *not_supported,
            role="null",
        ),
        KangCase(
            "fake_split_isg15",
            "False: ISG15 differs between fake donor groups (control cells only)",
            split,
            split_path,
            _de(
                "FALSE: ISG15 differs in CD14+ monocytes between fake groups A and B "
                "(Kang ctrl only)",
                MONO,
                ["ISG15"],
                "either",
                fake_cols,
                "A",
                "B",
            ),
            "NOT SUPPORTED",
            *not_supported,
            role="false",
        ),
        KangCase(
            "fake_split_lyz",
            "False: LYZ differs between fake donor groups (control cells only)",
            split,
            split_path,
            _de(
                "FALSE: LYZ differs in CD14+ monocytes between fake groups A and B "
                "(Kang ctrl only)",
                MONO,
                ["LYZ"],
                "either",
                fake_cols,
                "A",
                "B",
            ),
            "NOT SUPPORTED",
            *not_supported,
            role="false",
        ),
        KangCase(
            "fake_split_composition",
            "False: CD14+ monocyte proportion differs between fake donor groups",
            split,
            split_path,
            _spec(
                claim_type="composition_shift",
                name="FALSE: CD14+ monocyte proportion differs between fake groups A and B "
                "(Kang ctrl only)",
                columns=fake_cols,
                target=MONO,
                comparison={"group_a": "A", "group_b": "B", "direction": "different"},
            ),
            "NOT SUPPORTED",
            *not_supported,
            role="false",
        ),
        KangCase(
            "fake_shuffle_isg15",
            "False: ISG15 up in fake group A (condition shuffled at donor level)",
            shuffle,
            shuffle_path,
            _de(
                "FALSE: ISG15 up in CD14+ monocytes in fake group A vs B "
                "(Kang, donor-level shuffle)",
                MONO,
                ["ISG15"],
                "up_in_a",
                shuffle_cols,
                "A",
                "B",
            ),
            "NOT SUPPORTED",
            *not_supported,
            role="false",
        ),
    ]
    return cases


def load_kang(path: Path) -> AnnData:
    return ad.read_h5ad(path)
