"""Planted errors on PBMC3k.

Every case is built by a deterministic function of the raw data and a seed, states its claim,
and states *before running* what counts as "caught". The runner records the verdict and applies
the criterion mechanically; nothing here is adjusted after seeing results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from recluster.bench.pbmc3k import as_processed, expressed_genes
from recluster.config import GridConfig, PipelineParams
from recluster.data import detect_mito_prefix, extract_counts
from recluster.pipeline import Pipeline
from recluster.spec import ClaimSpec


@dataclass
class PlantedCase:
    key: str
    title: str
    how_built: str
    adata: AnnData
    spec: ClaimSpec
    audit: bool
    expected: str
    caught_rule: str
    caught: Any  # callable(result) -> bool
    extra: dict[str, Any] = field(default_factory=dict)
    role: str = "planted"  # or "control"


# ------------------------------------------------------------------ helpers


def _baseline_clusters(adata: AnnData, **overrides: Any) -> pd.Series:
    counts, _ = extract_counts(adata)
    pipe = Pipeline(counts, GridConfig(), mito_prefix=detect_mito_prefix(counts.var_names))
    params = PipelineParams(**{"batch_correction": False, **overrides})
    return pipe.run(params).clusters


def _all_holding_runs_have(result: Any, axis: str, value: Any) -> bool:
    holding = [r for r in result.runs if r.holds]
    return bool(holding) and all(r.params[axis] == value for r in holding)


def _top_breaker_axis(result: Any) -> str | None:
    return result.verdict.breakers[0]["axis"] if result.verdict.breakers else None


# ------------------------------------------------------------ 1. doublets


# names used in the case descriptions: (one cell, adjective for its markers)
_TYPE_WORDS = {
    "CD14 Mono": ("CD14 monocyte", "monocyte"),
    "B": ("B cell", "B-cell"),
    "NK": ("NK cell", "NK-cell"),
    "CD8 T": ("CD8 T cell", "CD8 T-cell"),
}


def _type_words(cell_type: str) -> tuple[str, str]:
    return _TYPE_WORDS.get(cell_type, (f"{cell_type} cell", cell_type))


def plant_doublets(
    adata: AnnData,
    reference: pd.Series,
    *,
    type_a: str = "CD14 Mono",
    type_b: str = "B",
    n: int = 150,
    seed: int = 0,
    markers: tuple[str, ...] = ("CD14", "LYZ", "MS4A1", "CD79A"),
    key: str = "doublets",
    name: str | None = None,
    held_out: bool = False,
) -> list[PlantedCase]:
    """Sum counts of random ``type_a`` + ``type_b`` pairs and append them as new "cells".

    A naive analysis sees a cluster carrying both types' markers (by default CD14+ LYZ+ MS4A1+
    CD79A+) and calls it a population. Each source cell is used at most once per type.
    ``held_out`` marks the cases pre-registered in PLAN.md "v0.6 plan" (H1, H2).
    """
    rng = np.random.default_rng(seed)
    a_cells = reference.index[reference == type_a]
    b_cells = reference.index[reference == type_b]
    for t, cells in ((type_a, a_cells), (type_b, b_cells)):
        if len(cells) < n:
            raise ValueError(f"{n} doublets need {n} '{t}' cells; the reference has {len(cells)}")
    pa = rng.choice(a_cells, size=n, replace=False)
    pb = rng.choice(b_cells, size=n, replace=False)
    x = sp.csr_matrix(adata.X)
    ia = adata.obs_names.get_indexer(pa)
    ib = adata.obs_names.get_indexer(pb)
    doublets = x[ia] + x[ib]
    new = AnnData(X=doublets, var=adata.var.copy())
    new.obs_names = [f"planted_doublet_{i}" for i in range(n)]
    out = AnnData(X=sp.vstack([x, doublets]).tocsr(), var=adata.var.copy())
    out.obs_names = list(adata.obs_names) + list(new.obs_names)
    out.obs["injected_doublet"] = np.r_[np.zeros(adata.n_obs, bool), np.ones(n, bool)]

    marker_text = " ".join(f"{g}+" for g in markers)
    prefix = "Planted doublets (held out)" if held_out else "Planted doublets"
    spec = ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": name or f"{prefix}: a {marker_text} population exists (PBMC3k)",
            "target": {"markers": list(markers), "min_cells": 20},
        }
    )
    (cell_a, adj_a), (cell_b, adj_b) = _type_words(type_a), _type_words(type_b)
    how_built = (
        f"{n} synthetic doublets made by summing the raw counts of a random {cell_a} "
        f"and a random {cell_b} (seed {seed}) were appended to PBMC3k. The claim asserts a "
        f"population carrying both {adj_a} and {adj_b} markers."
    )
    title = "Injected synthetic doublets forming a fake population"
    if held_out:
        title += f" (held out: {type_a} + {type_b})"
        how_built += (
            " Held-out case: its construction, claim and rule were fixed in PLAN.md before the "
            "doublet-population step was built, and it is run once with that step as adopted. "
            "If the doublets never form a population in any run the case is uninformative, and "
            "the rule scores it as missed."
        )
    extra: dict[str, Any] = {"n_injected": n}
    if held_out:
        extra.update({"held_out": True, "parent_types": [type_a, type_b], "seed": seed})
    return [
        PlantedCase(
            key=key,
            title=title,
            how_built=how_built,
            adata=out,
            spec=spec,
            audit=False,
            expected="not ROBUST; holds only when doublet removal is off",
            caught_rule=(
                "verdict is not ROBUST AND the claim held in at least one run AND every run in "
                "which it held had doublet_removal=off"
            ),
            caught=lambda r: (
                r.verdict.label != "ROBUST" and _all_holding_runs_have(r, "doublet_removal", False)
            ),
            extra=extra,
        )
    ]


def plant_held_out_doublets(adata: AnnData, reference: pd.Series) -> list[PlantedCase]:
    """The two held-out doublet cases pre-registered in PLAN.md "v0.6 plan".

    H1: 150 NK + CD14 monocyte doublets (seed 11). H2: 150 CD8 T + B doublets (seed 12).
    """
    return plant_doublets(
        adata,
        reference,
        type_a="NK",
        type_b="CD14 Mono",
        n=150,
        seed=11,
        markers=("GNLY", "NKG7", "CD14", "LYZ"),
        key="doublets_nk_mono",
        held_out=True,
    ) + plant_doublets(
        adata,
        reference,
        type_a="CD8 T",
        type_b="B",
        n=150,
        seed=12,
        markers=("CD8A", "CD8B", "MS4A1", "CD79A"),
        key="doublets_cd8_b",
        held_out=True,
    )


def scrublet_recall_on_injected(adata: AnnData) -> dict[str, float]:
    """Fraction of injected doublets Scrublet flags at the baseline QC setting."""
    counts, _ = extract_counts(adata)
    pipe = Pipeline(counts, GridConfig(), mito_prefix=detect_mito_prefix(counts.var_names))
    stage = pipe.qc_stage(GridConfig().baseline)
    obs = stage.adata.obs
    injected = obs["injected_doublet"].to_numpy()
    pred = obs["predicted_doublet"].to_numpy()
    return {
        "injected_passing_qc": int(injected.sum()),
        "recall": float(pred[injected].mean()) if injected.any() else float("nan"),
        "false_positive_rate": float(pred[~injected].mean()),
    }


# ----------------------------------------------------- 2. over-clustering


def plant_overclustering(
    adata: AnnData, reference: pd.Series, *, resolution: float = 3.0, seed: int = 0
) -> list[PlantedCase]:
    """Cluster at a high resolution and claim one shard of a homogeneous population.

    Rule for choosing the shard (fixed before running): among fine clusters that lie >= 90 %
    inside one reference population, take the parent population that is split into the most
    shards, and claim that parent's largest shard. The control is the fine cluster with the
    highest Jaccard to a whole reference population.
    """
    fine = _baseline_clusters(adata, resolution=resolution, seed=seed)
    ref = reference.loc[fine.index]
    table = pd.crosstab(fine, ref)
    purity = table.max(axis=1) / table.sum(axis=1)
    parent = table.idxmax(axis=1)
    pure = purity[purity >= 0.9].index
    splits = parent.loc[pure].value_counts()
    most_split = splits.idxmax()
    shards = [c for c in pure if parent[c] == most_split]
    spurious = max(shards, key=lambda c: int(table.loc[c].sum()))

    jaccards = {}
    for c in table.index:
        p = parent[c]
        inter = table.loc[c, p]
        union = table.loc[c].sum() + table[p].sum() - inter
        jaccards[c] = inter / union
    control = max(jaccards, key=lambda c: jaccards[c])

    labels = pd.Series("filtered", index=adata.obs_names)
    labels.loc[fine.index] = ("cluster_" + fine.astype(str)).to_numpy()
    spur_label, ctrl_label = f"cluster_{spurious}", f"cluster_{control}"
    processed = as_processed(adata, labels)
    processed.obs["reference"] = reference.loc[processed.obs_names].to_numpy()

    def spec(label: str, name: str) -> ClaimSpec:
        return ClaimSpec.model_validate(
            {
                "claim_type": "population_exists",
                "name": name,
                "columns": {"labels": "tool_label"},
                "target": {"label_value": label, "min_cells": 20},
            }
        )

    n_spur, n_ctrl = int(table.loc[spurious].sum()), int(table.loc[control].sum())
    built = (
        f"PBMC3k clustered with the baseline pipeline at Leiden resolution {resolution} gives "
        f"{len(table)} clusters. Reference population '{most_split}' is split into "
        f"{splits[most_split]} shards; the largest shard ({spur_label}, {n_spur} cells, "
        f"{purity[spurious]:.0%} '{most_split}') "
        f"is claimed as a population. Control: {ctrl_label} ({n_ctrl} cells), which coincides with "
        f"reference population '{parent[control]}' (Jaccard {jaccards[control]:.2f})."
    )
    return [
        PlantedCase(
            key="overclustering",
            title="Over-clustering that creates a spurious population",
            how_built=built,
            adata=processed,
            spec=spec(
                spur_label, f"Over-clustering: {spur_label} is a distinct population (PBMC3k)"
            ),
            audit=True,
            expected="not ROBUST",
            caught_rule="verdict is not ROBUST",
            caught=lambda r: r.verdict.label != "ROBUST",
            extra={"parent": most_split, "n_shards": int(splits[most_split]), "n_cells": n_spur},
        ),
        PlantedCase(
            key="overclustering_control",
            title="Control for over-clustering: a shard that is a whole population",
            how_built=f"Same object; {ctrl_label} ≈ reference '{parent[control]}'.",
            adata=processed,
            spec=spec(
                ctrl_label,
                f"Over-clustering control: {ctrl_label} is a distinct population (PBMC3k)",
            ),
            audit=True,
            expected="ROBUST",
            caught_rule="control passes if verdict is ROBUST",
            caught=lambda r: r.verdict.label == "ROBUST",
            extra={
                "parent": parent[control],
                "jaccard": float(jaccards[control]),
                "n_cells": n_ctrl,
            },
            role="control",
        ),
    ]


# ------------------------------------------------- 3. QC deletes a population


def plant_qc_deletion(
    adata: AnnData,
    reference: pd.Series,
    *,
    target: str = "NK",
    mito_low: float = 12.0,
    mito_high: float = 16.0,
    seed: int = 0,
) -> list[PlantedCase]:
    """Raise the mitochondrial fraction of one real population into the 12-16 % band.

    The baseline QC (max 10 % mito) then deletes it silently; only the 20 % setting keeps it.
    """
    rng = np.random.default_rng(seed)
    x = sp.lil_matrix(adata.X.astype(np.float32))
    prefix = detect_mito_prefix(adata.var_names) or "MT-"
    mt = np.flatnonzero(np.char.startswith(np.asarray(adata.var_names, dtype=str), prefix))
    cells = np.flatnonzero((reference.loc[adata.obs_names] == target).to_numpy())
    for i in cells:
        row = x.getrow(i).toarray().ravel()
        total, m = row.sum(), row[mt].sum()
        p = rng.uniform(mito_low, mito_high) / 100.0
        new_m = p / (1 - p) * (total - m)
        if m > 0:
            row[mt] = np.rint(row[mt] * new_m / m)
        else:
            row[mt[0]] = np.rint(new_m)
        x[i, mt] = row[mt]
    out = AnnData(X=x.tocsr(), obs=adata.obs.copy(), var=adata.var.copy())
    out.obs["reference"] = reference.loc[out.obs_names].to_numpy()

    spec = ClaimSpec.model_validate(
        {
            "claim_type": "population_exists",
            "name": f"QC deletion: {target} cells (GNLY+ NKG7+) exist (PBMC3k)",
            "target": {
                "markers": ["GNLY", "NKG7", "KLRD1"],
                "negative_markers": ["CD3D"],
                "min_cells": 20,
            },
        }
    )
    return [
        PlantedCase(
            key="qc_deletion",
            title="A QC threshold that silently deletes a real population",
            how_built=(
                f"Mitochondrial counts of the {len(cells)} reference {target} cells were rescaled "
                f"so each cell's mito fraction is uniform in [{mito_low:g}%, {mito_high:g}%] "
                f"(seed {seed}); all other "
                "counts untouched. The baseline max_mito_pct=10 removes every one of them."
            ),
            adata=out,
            spec=spec,
            audit=False,
            expected="not ROBUST; the failure table names max_mito_pct",
            caught_rule="top entry of the failure-rate table is the max_mito_pct axis",
            caught=lambda r: _top_breaker_axis(r) == "max_mito_pct",
            extra={"n_cells_modified": int(len(cells))},
        )
    ]


# ------------------------------------------- 4. confounded batch effect


def plant_batch_effect(
    adata: AnnData,
    reference: pd.Series,
    *,
    fraction: float = 0.4,
    n_genes: int = 300,
    strength: float = 2.5,
    seed: int = 0,
) -> list[PlantedCase]:
    """Multiply a random gene set in a random 40 % of cells; an uncorrected analysis splits
    populations by batch and calls one half a new subtype."""
    rng = np.random.default_rng(seed)
    x = sp.csr_matrix(adata.X.astype(np.float32)).tolil()
    batch2 = rng.random(adata.n_obs) < fraction
    up = expressed_genes(adata, 0.05, rng, n_genes)
    down = np.setdiff1d(expressed_genes(adata, 0.05, rng, 2 * n_genes), up)[:n_genes]
    rows = np.flatnonzero(batch2)
    sub = x[rows, :].toarray()
    sub[:, up] = np.rint(sub[:, up] * strength)
    sub[:, down] = np.rint(sub[:, down] / strength)
    x[rows, :] = sub
    out = AnnData(X=x.tocsr(), obs=adata.obs.copy(), var=adata.var.copy())
    out.obs["batch"] = np.where(batch2, "batch2", "batch1")
    out.obs["reference"] = reference.loc[out.obs_names].to_numpy()

    # the naive analysis: no batch correction
    clusters = _baseline_clusters(out, batch_correction=False)
    b = pd.Series(out.obs["batch"]).loc[clusters.index]
    ref = reference.loc[clusters.index]
    share2 = (b == "batch2").groupby(clusters.to_numpy()).mean()
    sizes = clusters.value_counts()
    artifact_clusters = share2[(share2 >= 0.9) & (sizes.loc[share2.index] >= 20)]
    if artifact_clusters.empty:
        raise RuntimeError(
            "batch planting did not produce a batch-pure cluster; strengthen the effect"
        )
    artifact = max(artifact_clusters.index, key=lambda c: int(sizes[c]))
    parent = ref[clusters == artifact].value_counts().idxmax()

    labels = pd.Series("filtered", index=out.obs_names)
    labels.loc[clusters.index] = ref.to_numpy()
    labels.loc[clusters.index[(clusters == artifact).to_numpy()]] = f"{parent} subtype"
    processed = as_processed(out, labels)

    def spec(label: str, name: str) -> ClaimSpec:
        return ClaimSpec.model_validate(
            {
                "claim_type": "population_exists",
                "name": name,
                "columns": {"labels": "tool_label", "batch": "batch"},
                "target": {"label_value": label, "min_cells": 20},
            }
        )

    n_art = int(sizes[artifact])
    built = (
        f"A random {fraction:.0%} of cells (seed {seed}) were assigned to 'batch2'; in those cells "
        f"{len(up)} random expressed genes were multiplied by {strength} and {len(down)} others "
        f"divided by {strength} (rounded). Clustering without correction yields cluster "
        f"{artifact} ({n_art} cells, {share2[artifact]:.0%} batch2, mostly '{parent}'), "
        f"labelled '{parent} subtype' and claimed as a "
        f"population. Control: the reference '{parent}' population itself."
    )
    return [
        PlantedCase(
            key="batch_effect",
            title="A fake batch effect confounded with a subset of cells",
            how_built=built,
            adata=processed,
            spec=spec(
                f"{parent} subtype",
                f"Batch artifact: '{parent} subtype' is a distinct population (PBMC3k)",
            ),
            audit=True,
            expected="not ROBUST; holds only when batch correction is off",
            caught_rule=(
                "verdict is not ROBUST AND every run in which the claim held had "
                "batch_correction=off"
            ),
            caught=lambda r: (
                r.verdict.label != "ROBUST" and _all_holding_runs_have(r, "batch_correction", False)
            ),
            extra={"artifact_cluster": str(artifact), "n_cells": n_art, "parent": parent},
        ),
        PlantedCase(
            key="batch_effect_control",
            title="Control for the batch effect: the real parent population",
            how_built=f"Same object; the reference '{parent}' cells as labelled.",
            adata=processed,
            spec=spec(parent, f"Batch control: '{parent}' is a distinct population (PBMC3k)"),
            audit=True,
            expected="ROBUST",
            caught_rule="control passes if verdict is ROBUST",
            caught=lambda r: r.verdict.label == "ROBUST",
            extra={"parent": parent},
            role="control",
        ),
    ]
