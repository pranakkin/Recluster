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
    # pairing is detected on the whole dataset; a run can still fall back to the unpaired
    # test when too few pairs have a pseudobulk in both groups, so the baseline's design is
    # reported separately
    paired = bool(data.get("paired_design"))
    baseline_design = data.get("baseline_design")

    p: list[str] = []
    p.append(
        f"Raw counts for {data['n_cells']:,} cells and {data['n_genes']:,} genes were taken "
        f"from {data['counts_source']} of the input file (SHA-256 {data['sha256'][:12]}…). "
        f"Each of {n_runs} pipeline runs (1 baseline, {grid_counts['one_at_a_time']} "
        f"one-at-a-time variations and {grid_counts['random']} seeded random combinations) "
        "re-derived clusters from these counts."
        + (
            f" {merged} grid point(s) that became identical to another after the QC cell-loss "
            "cap were dropped and are listed in the provenance."
            if (merged := grid_counts.get("merged_after_capping"))
            else ""
        )
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
        f"a threshold that alone would remove more than {grid.qc_max_cell_loss:.0%} of the cells "
        "was pulled back to the value that removes exactly that fraction, and the run records "
        "both values. Genes detected in fewer than 3 cells were dropped. Doublets were scored "
        f"with Scrublet (scanpy implementation{', per sample' if has_sample else ''}, seed "
        f"{grid.scrublet_seed}, exact neighbours); when its automatic threshold called fewer "
        f"than {grid.scrublet_fallback_trigger:.0%} of the expected doublets "
        f"({grid.scrublet_expected_rate_per_1000:.1%} per 1,000 cells) the top expected "
        "fraction by score was called instead. Predicted doublets were removed when the "
        f"doublet_removal axis was on (baseline: {'on' if b.doublet_removal else 'off'})."
        + _doublet_population_sentence(grid)
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
    pieces = []
    if tgt.markers:
        negative = (
            f" minus that of {', '.join(tgt.negative_markers)}" if tgt.negative_markers else ""
        )
        logic = (
            "the minimum of the per-marker z-scores, so every marker must be expressed"
            if tgt.marker_logic == "all"
            else f"the mean CP10k-log1p expression of {', '.join(tgt.markers)}{negative}"
        )
        where = (
            " within samples (a stratified within-sample AUROC: the Mann-Whitney U statistics "
            "of all samples pooled, so that the mix of other cell types a QC setting keeps "
            "cannot move it)"
            if has_sample
            else " against all other cells"
        )
        majority = (
            f" and more than {thr.marker_majority_min:.0%} of the markers were individually "
            f"enriched (per-marker AUROC at least {thr.marker_per_gene_auroc_min})"
            if len(tgt.markers) > 1 and tgt.marker_logic != "all"
            else ""
        )
        pieces.append(
            f"by marker enrichment: a per-cell score ({logic}, z-scored) was computed once from "
            f"the raw counts, and a cluster of at least {thr.min_cluster_size} cells qualified "
            f"when its AUROC for that score was at least {thr.marker_auroc_min}{where}{majority}"
        )
    if tgt.label_value is not None:
        pieces.append(
            f"by label overlap: a cluster of at least {thr.min_cluster_size} cells qualified when "
            f"at least {thr.label_precision_min:.0%} of its cells carry `{tgt.label_value}` in "
            f"`{spec.columns.labels}`"
        )
    joiner = " or " if len(pieces) == 2 and thr.candidate_rule == "either" else " and "
    p.append(
        "The target population was located in each run "
        + joiner.join(pieces)
        + "; the matched population is the union of qualifying clusters."
    )
    p.append(
        f"Overlap with the reference cell set ({reference_source}) was measured by the Jaccard "
        "index over barcodes that passed the run's QC; the fraction of reference cells removed "
        "by QC is reported alongside."
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
        if thr.composition_test == "logit_t":
            paired_test = "a paired t-test on the per-donor differences of logit-transformed "
            paired_test += "proportions"
            unpaired_test = "a Welch t-test on logit-transformed proportions"
        else:
            paired_test = "a Wilcoxon signed-rank test on per-donor differences"
            unpaired_test = "a Mann-Whitney U test on proportions"
        if paired:
            test = (
                f"{paired_test} when at least 2 donors have a sample in both groups, otherwise "
                f"{unpaired_test}"
            )
            test += _baseline_used(baseline_design, "the paired test", "the unpaired test")
        else:
            test = unpaired_test
        p.append(
            "For each run the proportion of cells in the population was computed per sample "
            f"(`{sample_key}`), and samples in `{cmp.group_a}` were compared with samples in "
            f"`{cmp.group_b}` by {test}, so that the biological replicate, not the cell, is "
            "the unit of inference. A run was only testable when the population had at least "
            f"{thr.min_cells_per_sample} cells in at least {thr.min_samples_present} samples of "
            f"each group. The claim held when p < {thr.alpha} and the direction was "
            f"`{cmp.direction}`. A cell-level chi-square test on pooled cells is reported only "
            "to show the effect of pseudoreplication."
        )
    else:
        cmp = spec.comparison
        assert cmp is not None
        if paired:
            design = (
                "~ pair + condition (blocked on the donor) when every pair has a pseudobulk in "
                "both groups, otherwise ~ condition"
            )
            design += _baseline_used(baseline_design, "~ pair + condition", "~ condition")
        else:
            design = "~ condition"
        p.append(
            "For each run, raw counts of cells in the population were summed per sample "
            f"(`{sample_key}`; samples with fewer than {thr.min_cells_per_sample} cells dropped) "
            f"and tested for `{cmp.group_a}` versus `{cmp.group_b}` with DESeq2 (pydeseq2; "
            f"genes with >= {thr.pseudobulk_min_total_counts} total counts). The design was "
            f"{design}. "
            f"A run was only testable with at least {thr.min_samples_present} pseudobulk "
            "replicates per group. "
            f"P-values for the {len(spec.genes)} claimed gene(s) were Benjamini-Hochberg "
            "adjusted over the claimed genes (genome-wide adjusted values are also shown). "
            f"The claim held when {'every' if spec.require == 'all' else 'any'} claimed gene "
            f"had adjusted p < {thr.alpha} with direction `{cmp.direction}`. A cell-level "
            "two-sided Wilcoxon rank-sum test on CP10k-log1p values is reported alongside; a "
            "gene significant only at the cell level is flagged as likely pseudoreplication "
            "(Squair et al., Nat Commun 2021)."
        )
    if spec.comparison is not None:
        # the permutation test enumerates every label assignment when there are few enough,
        # and falls back to seeded random draws otherwise; its unit is the sample unless samples
        # nest in donors, in which case the donor is the unit
        permutation = (
            "a permutation test (exact when at most 20,000 label assignments exist, otherwise "
            "5,000 seeded draws)"
        )
        units = (
            "at the level of the sample, or of the donor when samples are nested in donors that "
            "share one condition and one collection time"
        )
        if spec.columns.collection_time:
            p.append(
                "Time of day: each sample's collection time was converted to an angle, and the "
                f"association between time and condition was tested with {permutation} on the "
                f"distance between the groups' mean (cos, sin) vectors, {units}; the claim's "
                "replicate-level model was then refitted with cos(2πt/24) and sin(2πt/24) as "
                "covariates to see whether the effect survives."
            )
        else:
            p.append(
                "Time of day: with no collection times, the per-sample expression of 13 core "
                "clock genes (Takahashi 2017; Wucher et al. 2023) was summarised as the "
                "difference between the E-box arm and the activator arm of the clock and "
                f"compared between groups with {permutation}, {units}; this is a screen for a "
                "collection-time confound, not an estimate of circadian phase."
            )
    p.append(
        f"Verdicts: ROBUST if the claim held in at least {thr.robust_min_survival:.0%} of runs, "
        f"FRAGILE if in at least {thr.fragile_min_survival:.0%}, otherwise NOT SUPPORTED; a "
        "possible time confound lowers the verdict as described in the time-of-day section. "
        "All thresholds are listed below and are configurable."
    )
    return " ".join(p)


def _doublet_population_sentence(grid: GridConfig) -> str:
    """The cluster-level doublet removal that follows clustering in runs with doublet removal on.

    The minimum of 10 cells is the default of ``doublet_population_test``, which the pipeline uses.
    """
    if not grid.doublet_cluster_removal:
        return " Clusters were not tested as doublet populations (doublet_cluster_removal off)."
    return (
        " In those runs, after clustering, every cluster of at least 10 cells was compared with "
        "simulated doublets of two other cell-type groups (the other clusters grouped by the "
        "correlation of their mean expression; each doublet sums the raw counts of one random "
        "cell from each group, for the 3 pairs of groups that best explain the cluster's mean "
        "expression), and a cluster that a cross-validated logistic-regression classifier could "
        "not separate from them (out-of-fold AUROC, mean of 5 repetitions, below "
        f"{grid.doublet_cluster_auroc_max:.2f} for some pair) was removed as a doublet "
        "population."
    )


def _opts(values: list[Any], baseline: Any) -> str:
    shown = [f"{v}" + ("*" if v == baseline else "") for v in values]
    return "{" + ", ".join(shown) + "}"


def _baseline_used(design: Any, paired: str, unpaired: str) -> str:
    """Clause naming the design the baseline run used, or nothing when it is not recorded."""
    if design == "paired":
        return f"; the baseline run used {paired}"
    if design == "unpaired":
        return f"; the baseline run used {unpaired}"
    return ""
