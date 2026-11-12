"""Kang et al. 2018 (Nat Biotechnol 36:89): PBMCs from 8 lupus patients, IFN-beta vs control.

Source: the processed object distributed with pertpy (``kang_2018.h5ad`` on
exampledata.scverse.org), raw integer counts in .X, columns ``label`` (ctrl/stim),
``replicate`` (8 donors), ``cell_type``.
This object contains singlets only; the demuxlet multiplet calls from GEO are *not* in it.
"""

from __future__ import annotations

import hashlib
import math
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
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

# Planted time-of-day signal. These values are used only to PLANT a signal; nothing here is
# estimated from Kang. Peak hours are rounded values consistent with the human blood and tissue
# literature (Ruben et al., Sci Transl Med 2018, human circadian atlas): the E-box arm
# PER1/PER3/NR1D1/NR1D2/DBP/TEF/CIART peaks in the morning, PER2 around midday, CRY1 is delayed
# to the evening (it peaks later than the other E-box targets: Ukai-Tadenuma et al., Cell 2011),
# and the activator arm ARNTL/NPAS2/CLOCK peaks in the evening/night. Half-amplitudes of
# 0.2-0.7 log2 are 1.3-2.6-fold peak to trough, within the range reported for core clock genes
# in human tissues.
#
# Day/night direction was checked against the Wucher et al. 2023 S8 table shipped as
# recluster/resources/daynight_genes.csv (``time_signal_day_night``: a morning or midday peak,
# <= 14 h, should be day-higher and an evening peak, >= 18 h, night-higher). 11 of the 14 genes
# agree. Three do not: CRY1 is day-higher in 20 of 20 GTEx tissues, NR1D1 night-higher in 22 of
# 22, and NR1D2 has no consistent direction (9 day vs 11 night, binomial P 0.82). The S8 table
# splits post-mortem donors by time of death into day and night and does not resolve a peak
# hour, so the peaks are kept as documented above and the disagreement is recorded in
# ``TIME_SIGNAL_S8_DISAGREE`` rather than tuned away; the test suite checks that the list is
# exactly what the shipped table gives.
#
# SQSTM1 is the claimed gene: higher by night in 6 of 6 GTEx tissues (S8), well expressed in
# monocytes, with a donor-to-donor log2 CPM range of 0.4 across Kang's 8 control donors within
# CD14+ monocytes (0.6 over all control cells). It replaced RPS26 (also on the S8 list, day-higher
# in 16 of 16 tissues), whose range across the control donors is about 2.6 log2 (1.7 within
# CD14+ monocytes) and masked a 0.4 log2 planted shift. That choice is a documented design
# decision about which gene carries the planted signal, not a tuning of the signal.
TIME_SIGNAL: dict[str, tuple[float, float]] = {  # gene: (peak hour, half-amplitude log2)
    "PER1": (7, 0.5),
    "PER2": (14, 0.4),
    "PER3": (7, 0.7),
    "CRY1": (20, 0.25),  # delayed E-box target; S8 direction disagrees, see above
    "CRY2": (10, 0.3),
    "NR1D1": (7, 0.7),  # S8 direction disagrees, see above
    "NR1D2": (8, 0.6),  # S8 has no consistent direction, see above
    "DBP": (8, 0.7),
    "TEF": (9, 0.6),
    "CIART": (8, 0.6),
    "ARNTL": (20, 0.4),
    "NPAS2": (20, 0.4),
    "CLOCK": (20, 0.2),
    "SQSTM1": (20, 0.6),  # 2.3-fold peak to trough: a strong but documented rhythmic amplitude
}
# Genes whose TIME_SIGNAL peak does not match the day/night direction in the S8 table.
TIME_SIGNAL_S8_DISAGREE = ("CRY1", "NR1D1", "NR1D2")
TIME_CLAIM_GENE = "SQSTM1"
TIME_CLAIM_GROUP_A = "B"  # the afternoon group, where an evening-peaking gene is higher
TIME_CLAIM_GROUP_B = "A"
MORNING_TIMES = [8.0, 9.0, 10.5, 12.0]
AFTERNOON_TIMES = [13.0, 15.0, 16.5, 18.0]


def time_signal_day_night(table: pd.DataFrame | None) -> dict[str, dict[str, Any]]:
    """Compare each TIME_SIGNAL peak with the day/night direction in the Wucher 2023 S8 table.

    A morning or midday peak (<= 14 h) should be day-higher (more day than night tissues) and an
    evening peak (>= 18 h) night-higher. Returns one row per gene present in the table with the
    expected and observed direction and whether they agree; genes absent from the table are
    left out. The rule exists so that the provenance comment above can be checked mechanically.
    """
    if table is None or table.empty:
        return {}
    rows = table.drop_duplicates("gene").set_index("gene")
    out: dict[str, dict[str, Any]] = {}
    for gene, (peak, _amp) in TIME_SIGNAL.items():
        if gene not in rows.index:
            continue
        day, night = int(rows.loc[gene, "day_tissues"]), int(rows.loc[gene, "night_tissues"])
        expected = "day" if peak <= 14 else "night" if peak >= 18 else "unspecified"
        observed = "day" if day > night else "night" if night > day else "tie"
        out[gene] = {
            "peak_hour": float(peak),
            "expected": expected,
            "observed": observed,
            "day_tissues": day,
            "night_tissues": night,
            "agrees": expected == observed,
        }
    return out


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
    role: str = "true"  # true | null | false | planted_time


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


def plant_time_signal(adata: AnnData, *, seed: int = 0) -> AnnData:
    """All cells; donors split into fake groups A/B with group A collected in the morning.

    Each donor keeps one collection time for both of its samples (the ctrl and stim halves were
    processed together). For every gene in ``TIME_SIGNAL`` the counts of all cells in a sample
    are multiplied by ``2 ** (amplitude * cos(2 * pi * (t - peak) / 24))`` and Poisson-resampled
    so that the injected effect is a smooth function of time, not a group label.
    """
    rng = np.random.default_rng(seed)
    out = adata.copy()
    donors = sorted(out.obs["replicate"].astype(str).unique())
    rng.shuffle(donors)
    half = len(donors) // 2
    group = {d: ("A" if i < half else "B") for i, d in enumerate(donors)}
    times = {}
    for i, d in enumerate(donors):
        times[d] = MORNING_TIMES[i] if i < half else AFTERNOON_TIMES[i - half]
    rep = out.obs["replicate"].astype(str)
    out.obs["fake_group"] = pd.Categorical(rep.map(group))
    out.obs["collection_time"] = rep.map(times).astype(float)
    x = sp.lil_matrix(out.X.astype(np.float32))
    factors: dict[str, dict[str, float]] = {}
    for gene, (peak, amp) in TIME_SIGNAL.items():
        if gene not in out.var_names:
            continue
        j = int(out.var_names.get_loc(gene))
        col = np.asarray(out.X[:, j].todense()).ravel() if sp.issparse(out.X) else out.X[:, j]
        new = col.copy()
        factors[gene] = {}
        for d, t in times.items():
            f = 2.0 ** (amp * math.cos(2 * math.pi * (t - peak) / 24.0))
            factors[gene][d] = float(f)
            rows = np.flatnonzero((rep == d).to_numpy())
            new[rows] = rng.poisson(col[rows] * f)
        x[:, j] = new[:, None]
    out.X = x.tocsr()
    out.uns["planted_time"] = {
        "groups": group,
        "collection_time_by_donor": times,
        "factors_by_gene_and_donor": factors,
        "genes": sorted(factors),
    }
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


def _held_unadjusted(tc: dict[str, Any] | None) -> bool | None:
    """The adjustment's held_unadjusted flag from the full or the compacted time-check dict."""
    if not tc:
        return None
    if "held_unadjusted" in tc:
        return tc["held_unadjusted"]
    return (tc.get("adjustment") or {}).get("held_unadjusted")


def kang_cases(
    real_path: Path,
    split_path: Path,
    shuffle_path: Path,
    split_groups: dict[str, list[str]],
    time_path: Path | None = None,
    time_groups: dict[str, str] | None = None,
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
    if time_path is not None:
        tg = time_groups or {}
        a_donors = sorted(d for d, g in tg.items() if g == "A")
        b_donors = sorted(d for d, g in tg.items() if g == "B")
        planted = (
            "All cells. The 8 donors were split (seed 0) into fake groups A = "
            + ", ".join(a_donors)
            + " (collected 08:00-12:00) and B = "
            + ", ".join(b_donors)
            + " (13:00-18:00), 4 vs 4. Each donor is one sample: its ctrl and stim cells are "
            + "pooled (they were processed together and share one collection time), so the "
            + "analysis has 8 units, not 16. Counts of 13 core clock genes and of SQSTM1 (higher "
            + "by night in 6 of 6 GTEx tissues, Wucher 2023) were multiplied per donor by "
            + "2^(a cos(2π(t - peak)/24)) with the documented peaks and half-amplitudes "
            + "(0.2-0.7 log2) and Poisson-resampled, so the only difference between A and B is a "
            + "time-of-day effect. IFN-beta status is balanced within every donor and cannot "
            + "confound the fake groups."
        )
        time_cols_b = {"condition": "fake_group", "sample": ["replicate"], "labels": "cell_type"}
        time_cols_a = dict(time_cols_b, collection_time="collection_time")
        cases.append(
            KangCase(
                "time_mode_b",
                "Planted time-of-day effect claimed as a condition effect (times unknown)",
                planted,
                time_path,
                _de(
                    f"PLANTED TIME, times unknown: {TIME_CLAIM_GENE} up in CD14+ monocytes, "
                    "afternoon group B vs morning group A (Kang)",
                    MONO,
                    [TIME_CLAIM_GENE],
                    "up_in_a",
                    time_cols_b,
                    TIME_CLAIM_GROUP_A,
                    TIME_CLAIM_GROUP_B,
                ),
                "label before the time check ROBUST; time check POSSIBLE TIME CONFOUND; "
                "final verdict not ROBUST",
                "label_before_time_check is ROBUST (the replicate-aware analysis called the "
                "planted effect) AND the time-check label is POSSIBLE TIME CONFOUND AND the "
                "final verdict is not ROBUST",
                lambda r: (
                    r.verdict.label_before_time_check == "ROBUST"
                    and (r.time_check or {}).get("label") == "POSSIBLE TIME CONFOUND"
                    and r.verdict.label != "ROBUST"
                ),
                role="planted_time",
            )
        )
        cases.append(
            KangCase(
                "time_mode_a",
                "Same planted effect with collection times given",
                planted + " Collection times are passed to Recluster.",
                time_path,
                _de(
                    f"PLANTED TIME, times given: {TIME_CLAIM_GENE} up in CD14+ monocytes, "
                    "afternoon group B vs morning group A (Kang)",
                    MONO,
                    [TIME_CLAIM_GENE],
                    "up_in_a",
                    time_cols_a,
                    TIME_CLAIM_GROUP_A,
                    TIME_CLAIM_GROUP_B,
                ),
                "claim held before adjustment; effect removed by the time adjustment; "
                "verdict NOT SUPPORTED",
                "the time check's adjustment has held_unadjusted True (the claim held before "
                "adjusting for collection time) AND effect_removed True (it is not significant "
                "after adjusting) AND the verdict is NOT SUPPORTED",
                lambda r: (
                    _held_unadjusted(r.time_check) is True
                    and (r.time_check or {}).get("effect_removed") is True
                    and r.verdict.label == "NOT SUPPORTED"
                ),
                role="planted_time",
            )
        )
    return cases


def load_kang(path: Path) -> AnnData:
    return ad.read_h5ad(path)
