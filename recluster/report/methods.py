"""Deterministic methods paragraph. A template filled with the actual settings; no LLM."""

from __future__ import annotations

from typing import Any

from recluster.config import GridConfig, Thresholds
from recluster.spec import ClaimSpec


def methods_paragraph(
    spec: ClaimSpec,
    grid: GridConfig,
    thr: Thresholds,
    data: dict[str, Any],
    n_runs: int,
    grid_counts: dict[str, int],
    reference_source: str,
) -> str:
    ax = grid.axes
    b = grid.baseline
    has_batch = spec.columns.batch is not None
    has_sample = spec.columns.sample is not None
    sample_key = "`, `".join(spec.columns.sample_key or [])

    p: list[str] = []
    p.append(
        f"Raw counts for {data['n_cells']:,} cells and {data['n_genes']:,} genes were taken "
        f"from {data['counts_source']} of the input file (SHA-256 {data['sha256'][:12]}…). "
        f"Each of {n_runs} pipeline runs (1 baseline, {grid_counts['one_at_a_time']} "
        f"one-at-a-time variations and {grid_counts['random']} seeded random combinations) "
        "re-derived clusters from these counts."
    )
    mito = (
        f"and had at most {_opts(ax.max_mito_pct, b.max_mito_pct)}% mitochondrial counts; "
        if data.get("mito_prefix") is not None
        else (
            "(the input has no mitochondrial genes, so the mitochondrial-fraction axis "
            "was not varied); "
        )
    )
    p.append(
        f"Cells were kept if they expressed at least {_opts(ax.min_genes, b.min_genes)} genes "
        f"{mito}"
        "genes detected in fewer than 3 cells were dropped. Doublets were scored with Scrublet "
        f"(scanpy implementation{', per sample' if has_sample else ''}, seed "
        f"{grid.scrublet_seed}) and removed when the doublet_removal axis was on "
        f"(baseline: {'on' if b.doublet_removal else 'off'})."
    )
    harmony = (
        f", optionally corrected for the batch column `{spec.columns.batch}` with Harmony "
        f"(baseline: {'on' if b.batch_correction else 'off'})"
        if has_batch
        else ""
    )
    p.append(
        "Expression was normalised either by total-count scaling to 10,000 counts per cell "
        f"followed by log1p and z-scaling of {_opts(ax.n_hvgs, b.n_hvgs)} highly variable genes "
        "(Seurat flavour), or by analytic Pearson residuals on the same number of genes "
        f"(baseline: {b.normalization}). {_opts(ax.n_pcs, b.n_pcs)} principal components "
        f"(ARPACK, seeded) were computed{harmony}; a {grid.n_neighbors}-nearest-neighbour graph "
        "(exact search) was clustered with Leiden (igraph flavour, resolution "
        f"{_opts(ax.resolution, b.resolution)}, seeds {_opts(ax.seed, b.seed)})."
    )

    tgt = spec.target
    if tgt.markers:
        negative = (
            f" minus that of {', '.join(tgt.negative_markers)}" if tgt.negative_markers else ""
        )
        p.append(
            "The target population was located in each run by marker enrichment: a per-cell "
            f"score (mean CP10k-log1p expression of {', '.join(tgt.markers)}{negative}, "
            "z-scored) was computed once from the raw counts, and clusters with at least "
            f"{thr.min_cluster_size} cells whose AUROC for that score versus all other cells was "
            f"at least {thr.marker_auroc_min} were taken as the population."
        )
    else:
        p.append(
            "The target population was located in each run as the union of clusters (at least "
            f"{thr.min_cluster_size} cells) in which at least {thr.label_precision_min:.0%} of "
            f"cells carry the label `{tgt.label_value}` in `{spec.columns.labels}`."
        )
    p.append(
        f"Overlap with the reference cell set ({reference_source}) was measured by the Jaccard "
        "index over barcodes, counting cells removed by QC as absent."
    )

    if spec.claim_type == "population_exists":
        replicate = (
            f", and at least {thr.min_cells_per_sample} cells in each of at least "
            f"{thr.min_samples_present} samples."
            if has_sample
            else "."
        )
        marker = f", marker AUROC >= {thr.marker_auroc_min}" if tgt.markers else ""
        p.append(
            f"The claim held in a run when the matched population had at least {tgt.min_cells} "
            f"cells, Jaccard >= {thr.jaccard_min}{marker}{replicate}"
        )
    elif spec.claim_type == "composition_shift":
        cmp = spec.comparison
        assert cmp is not None
        test = (
            "Welch t-test on logit-transformed proportions"
            if thr.composition_test == "logit_t"
            else "Mann-Whitney U test on proportions"
        )
        p.append(
            "For each run the proportion of cells in the population was computed per sample "
            f"(`{sample_key}`), and samples in `{cmp.group_a}` were compared with samples in "
            f"`{cmp.group_b}` by a {test}, so that the biological replicate, not the cell, is "
            f"the unit of inference. The claim held when p < {thr.alpha} and the direction was "
            f"`{cmp.direction}`. A cell-level chi-square test on pooled cells is reported only "
            "to show the effect of pseudoreplication."
        )
    else:
        cmp = spec.comparison
        assert cmp is not None
        p.append(
            "For each run, raw counts of cells in the population were summed per sample "
            f"(`{sample_key}`; samples with fewer than {thr.min_cells_per_sample} cells dropped) "
            f"and tested for `{cmp.group_a}` versus `{cmp.group_b}` with DESeq2 (pydeseq2, "
            f"design ~condition, genes with >= {thr.pseudobulk_min_total_counts} total counts). "
            f"P-values for the {len(spec.genes)} claimed gene(s) were Benjamini-Hochberg "
            "adjusted over the claimed genes (genome-wide adjusted values are also shown). "
            f"The claim held when {'every' if spec.require == 'all' else 'any'} claimed gene "
            f"had adjusted p < {thr.alpha} with direction `{cmp.direction}`. A cell-level "
            "two-sided Wilcoxon rank-sum test on CP10k-log1p values is reported alongside; a "
            "gene significant only at the cell level is flagged as likely pseudoreplication "
            "(Squair et al., Nat Commun 2021)."
        )
    p.append(
        f"Verdicts: ROBUST if the claim held in at least {thr.robust_min_survival:.0%} of runs, "
        f"FRAGILE if in at least {thr.fragile_min_survival:.0%}, otherwise NOT SUPPORTED. "
        "All thresholds are listed below and are configurable."
    )
    return " ".join(p)


def _opts(values: list[Any], baseline: Any) -> str:
    shown = [f"{v}" + ("*" if v == baseline else "") for v in values]
    return "{" + ", ".join(shown) + "}"
