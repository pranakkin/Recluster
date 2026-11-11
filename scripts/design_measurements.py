"""Recompute the measurements cited for the design decisions of v0.2 to v0.5.

bench/RESULTS.md (section "What changed between the runs") and docs/HOW_IT_WORKS.md (section 8)
quote numbers that were measured while the design was being settled. This script recomputes
each of them from the files in data/ with the package's own functions and writes
bench/design_measurements.json (every number, unrounded) and bench/design_measurements.md (a
table of cited against measured values). Everything is seeded; a second run gives the same
numbers.

Inputs: data/kang_2018.h5ad, data/geo/ (GEO GSE96583 batch 2), data/planted/pbmc3k_doublets.h5ad
and, for item 3, bench/out/kang/ from `recluster bench --part kang`.

Usage:
    python scripts/design_measurements.py [--data data] [--bench bench] [--only 1,2,...]
        [--cpus 3]

The cluster-level doublet mixture analysis (RESULTS item 9) is not recomputed here.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

import recluster
from recluster.bench import kang as kang_module
from recluster.bench.doublets import _agreement, load_batch2
from recluster.bench.kang import load_kang
from recluster.config import GridConfig, PipelineParams, Thresholds
from recluster.data import detect_mito_prefix, extract_counts, load_adata
from recluster.evaluators.composition import evaluate_composition
from recluster.evaluators.de import pseudobulk_de
from recluster.evaluators.pairing import pair_key
from recluster.evaluators.pseudobulk import pseudobulk_counts
from recluster.matching import (
    auroc_by_group,
    auroc_two_sets,
    lognorm_matrix,
    marker_score,
    within_sample_auroc,
)
from recluster.pipeline import SAMPLE_COL, Pipeline
from recluster.spec import Comparison, make_sample_key
from recluster.time_confound import exchangeable_units, panel_check

ROOT = Path(__file__).resolve().parents[1]
CD4T_CASE = "kang/isg15-up-in-cd4-t-cells-stim-vs-ctrl-kang-2018/provenance.json"
MONO_LABEL = "CD14+ Monocytes"
MONO_MARKERS = ["CD14", "LYZ", "S100A8"]  # as in recluster.bench.kang.MONO
SAMPLE_KEY = ["replicate", "label"]  # donor x condition, as in recluster.bench.kang.COLS
CONDITION = "label"

Row = dict[str, Any]


def say(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def row(
    item: str,
    quantity: str,
    cited: str,
    measured: str,
    how: str,
    agrees: bool | None,
) -> Row:
    """One line of the output table. ``agrees`` is None when the comparison is not numeric."""
    return {
        "item": item,
        "quantity": quantity,
        "cited": cited,
        "measured": measured,
        "agrees": agrees,
        "how": how,
    }


def same_at(measured: float, cited: float, digits: int) -> bool:
    """Does ``measured`` round to the cited value at the cited precision?"""
    return round(measured, digits) == round(cited, digits)


class Kang:
    """The Kang 2018 counts and the keys every Kang item uses, loaded once."""

    def __init__(self, data_dir: Path) -> None:
        self.adata = load_kang(data_dir / "kang_2018.h5ad")
        self.counts, _ = extract_counts(self.adata)
        obs = self.counts.obs
        self.sample = make_sample_key(obs, SAMPLE_KEY)  # donor x condition, 16 samples
        self.condition = obs["label"].astype(str)
        self.labels = obs["cell_type"].astype(str)
        self.donor = obs["replicate"].astype(str)


# --------------------------------------------------------------------------- item 1


def item_qc_cap(kang: Kang, cpus: int) -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 2: what min_genes=500 does to Kang, with and without the cell-loss cap.

    The stimulated-monocyte cluster of a run is the cluster holding the most cells labelled
    CD14+ Monocytes in the stim condition. Its whole-dataset marker AUROC is the statistic the
    first benchmark used (``auroc_by_group`` on the CD14/LYZ/S100A8 score over the run's cells).
    The first run predates the Scrublet fallback, so the two runs are repeated with the
    automatic threshold alone (fallback trigger 0) to show how much that difference matters.
    """
    counts = kang.counts
    n_genes = np.asarray((counts.X > 0).sum(axis=1)).ravel()
    removed = float((n_genes < 500).mean())
    mito = detect_mito_prefix(counts.var_names)
    capped_pipe = Pipeline(counts, GridConfig(), sample=kang.sample, mito_prefix=mito)
    capped_value, _, _ = capped_pipe.effective_qc(PipelineParams(min_genes=500))

    lognorm = lognorm_matrix(counts, MONO_MARKERS)
    score = marker_score(lognorm, MONO_MARKERS, [])
    stim_mono = (kang.labels == MONO_LABEL) & (kang.condition == "stim")

    variants = {
        "current pipeline": GridConfig(qc_max_cell_loss=1.0, n_cpus=cpus),
        "automatic Scrublet threshold only": GridConfig(
            qc_max_cell_loss=1.0, scrublet_fallback_trigger=0.0, n_cpus=cpus
        ),
    }
    runs: dict[str, dict[str, Any]] = {}
    for name, cfg in variants.items():
        pipe = Pipeline(counts, cfg, sample=kang.sample, mito_prefix=mito)
        runs[name] = {}
        for min_genes in (200, 500):
            say(f"  item 1: {name}, min_genes={min_genes}")
            rr = pipe.run(PipelineParams(min_genes=min_genes, batch_correction=False))
            cl = rr.clusters
            target = str(cl[stim_mono.loc[cl.index].to_numpy()].value_counts().idxmax())
            in_target = (cl == target).to_numpy()
            s = score.loc[cl.index].to_numpy()
            overall = auroc_by_group(s, cl.to_numpy())[target]
            within = within_sample_auroc(s, cl.to_numpy(), kang.sample.loc[cl.index].to_numpy())
            kept_labels = kang.labels.loc[cl.index]
            runs[name][str(min_genes)] = {
                "n_cells_kept": int(len(cl)),
                "n_removed_qc": rr.n_removed_qc,
                "n_removed_doublet": rr.n_removed_doublet,
                "n_clusters": rr.n_clusters,
                "stim_mono_cluster": target,
                "cluster_size": int(in_target.sum()),
                "cluster_share_stim_mono": float(
                    stim_mono.loc[cl.index].to_numpy()[in_target].mean()
                ),
                "whole_dataset_auroc": overall,
                "within_sample_auroc": within[target],
                "mono_share_of_kept_cells": float((kept_labels == MONO_LABEL).mean()),
                "lymphocyte_share_of_kept_cells": float(
                    kept_labels.isin(["CD4 T cells", "CD8 T cells", "B cells", "NK cells"]).mean()
                ),
            }
    cur = runs["current pipeline"]
    old = runs["automatic Scrublet threshold only"]
    how_auroc = (
        "Pipeline.run at min_genes 200 and 500 with GridConfig(qc_max_cell_loss=1.0); cluster "
        "with the most stim cells labelled CD14+ Monocytes; auroc_by_group of the CD14/LYZ/S100A8 "
        "marker_score over the run's cells"
    )
    rows = [
        row(
            "1 (RESULTS 2)",
            "Stim-monocyte cluster whole-dataset AUROC, baseline -> min_genes=500 uncapped",
            "0.83 -> 0.79",
            f"{cur['200']['whole_dataset_auroc']:.3f} -> {cur['500']['whole_dataset_auroc']:.3f}"
            f" (automatic Scrublet threshold only: {old['200']['whole_dataset_auroc']:.3f} -> "
            f"{old['500']['whole_dataset_auroc']:.3f})",
            how_auroc,
            same_at(cur["200"]["whole_dataset_auroc"], 0.83, 2)
            and same_at(cur["500"]["whole_dataset_auroc"], 0.79, 2),
        ),
        row(
            "1 (RESULTS 2)",
            "Same cluster, stratified within-sample AUROC (the current statistic)",
            "not cited",
            f"{cur['200']['within_sample_auroc']:.3f} -> {cur['500']['within_sample_auroc']:.3f}",
            "within_sample_auroc on the same runs, sample = replicate x label",
            None,
        ),
        row(
            "1 (RESULTS 2)",
            "Fraction of Kang cells min_genes=500 removes",
            "44 %",
            f"{100 * removed:.1f} %",
            "share of cells with fewer than 500 detected genes in the raw counts",
            same_at(100 * removed, 44, 0),
        ),
        row(
            "1 (RESULTS 2)",
            "min_genes=500 after the default 25 % cap",
            "441",
            str(capped_value),
            "Pipeline.effective_qc(PipelineParams(min_genes=500)) with GridConfig()",
            capped_value == 441,
        ),
    ]
    details = {
        "fraction_removed_by_min_genes_500": removed,
        "median_genes_per_cell": float(np.median(n_genes)),
        "capped_min_genes": capped_value,
        "runs": runs,
    }
    return rows, details


# --------------------------------------------------------------------------- item 2


def simulate_auroc_statistics(
    true_auroc: float,
    *,
    n_samples: int = 16,
    n_cluster: int = 12,
    n_other: int = 300,
    n_rep: int = 200,
    seed: int = 0,
) -> dict[str, Any]:
    """Stratified within-sample AUROC against the minimum of the per-sample AUROCs.

    Scores are binormal with equal variance: other cells N(0, 1), cluster cells N(d, 1) with
    d = sqrt(2) * Phi^-1(AUROC), so the population AUROC of every sample is ``true_auroc``.
    """
    rng = np.random.default_rng(seed)
    shift = math.sqrt(2.0) * float(stats.norm.ppf(true_auroc))
    per_sample = n_cluster + n_other
    samples = np.repeat(np.arange(n_samples), per_sample).astype(str)
    is_cluster = np.tile(np.r_[np.ones(n_cluster, bool), np.zeros(n_other, bool)], n_samples)
    groups = np.where(is_cluster, "cluster", "other")
    strat = np.empty(n_rep)
    minimum = np.empty(n_rep)
    for r in range(n_rep):
        score = rng.standard_normal(len(samples)) + shift * is_cluster
        strat[r] = within_sample_auroc(score, groups, samples)["cluster"]
        minimum[r] = min(
            auroc_two_sets(score[samples == s], is_cluster[samples == s])
            for s in np.unique(samples)
        )
    return {
        "true_auroc": true_auroc,
        "shift_d": shift,
        "seed": seed,
        "n_rep": n_rep,
        "stratified_mean": float(strat.mean()),
        "stratified_sd": float(strat.std(ddof=1)),
        "stratified_reject_below_0.80": float((strat < 0.80).mean()),
        "minimum_mean": float(minimum.mean()),
        "minimum_reject_below_0.80": float((minimum < 0.80).mean()),
    }


def item_simulation() -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 2: why the minimum over samples was not used."""
    sims = {
        "0.85": simulate_auroc_statistics(0.85, seed=0),
        "0.90": simulate_auroc_statistics(0.90, seed=1),
        "0.90, 2000 replicates": simulate_auroc_statistics(0.90, n_rep=2000, seed=2),
    }
    a, b, long_run = sims["0.85"], sims["0.90"], sims["0.90, 2000 replicates"]
    how = (
        "16 samples x (12 cluster + 300 other cells), binormal scores with equal variance, "
        "200 replicates (seed 0 at 0.85, seed 1 at 0.90); rejection = statistic < 0.80"
    )
    rows = [
        row(
            "2 (RESULTS 2)",
            "Rejection rate at true AUROC 0.85: minimum over samples vs stratified",
            "95 % vs 0 %",
            f"{100 * a['minimum_reject_below_0.80']:.1f} % vs "
            f"{100 * a['stratified_reject_below_0.80']:.1f} %",
            how,
            same_at(100 * a["minimum_reject_below_0.80"], 95, 0)
            and a["stratified_reject_below_0.80"] == 0,
        ),
        row(
            "2 (RESULTS 2)",
            "Mean of the stratified statistic at true AUROC 0.85",
            "0.851",
            f"{a['stratified_mean']:.3f}",
            "within_sample_auroc, mean over the 200 replicates",
            same_at(a["stratified_mean"], 0.851, 3),
        ),
        row(
            "2 (RESULTS 2)",
            "Rejection rate of the minimum at true AUROC 0.90, 200 replicates",
            "not cited (200 replicates vary; see the next row)",
            f"{100 * b['minimum_reject_below_0.80']:.1f} % (stratified: "
            f"{100 * b['stratified_reject_below_0.80']:.1f} %)",
            "as above, minimum of auroc_two_sets over the 16 samples",
            None,
        ),
        row(
            "2 (RESULTS 2)",
            "Rejection rate of the minimum at true AUROC 0.90, 2,000 replicates",
            "about 32 %",
            f"{100 * long_run['minimum_reject_below_0.80']:.1f} %",
            "as above, 2,000 replicates, seed 2",
            same_at(100 * long_run["minimum_reject_below_0.80"], 32, 0),
        ),
    ]
    return rows, sims


# --------------------------------------------------------------------------- item 3


def item_cd4t_jaccard(bench_dir: Path) -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 4: the CD4 T population under the labels-or-markers rule, from the bench."""
    path = bench_dir / "out" / CD4T_CASE
    if not path.exists():
        msg = f"{path} not found; run `recluster bench --part kang` first"
        say(f"  item 3: {msg}")
        rows = [
            row(
                "3 (RESULTS 4)",
                "CD4 T matched in every Kang run, Jaccard range",
                "every run, 0.90-0.97",
                f"not available: {msg}",
                "read from the benchmark's provenance.json",
                None,
            )
        ]
        return rows, {"available": False, "path": str(path)}
    prov = json.loads(path.read_text())
    runs = prov["runs"]
    found = [r for r in runs if r["match"]["found"]]
    jac = [float(r["match"]["jaccard"]) for r in found if r["match"]["jaccard"] is not None]
    rule = prov["thresholds"]["candidate_rule"]
    by = sorted({v for r in found for v in r["match"]["candidate_by"].values()})
    rows = [
        row(
            "3 (RESULTS 4)",
            "CD4 T matched in every Kang grid run, Jaccard range",
            "every run, 0.90-0.97",
            f"{len(found)} of {len(runs)} runs, {min(jac):.2f}-{max(jac):.2f}",
            f"bench/out/{CD4T_CASE}: runs[].match (candidate_rule={rule}; "
            f"candidates qualified by {', '.join(by)})",
            len(found) == len(runs) and same_at(min(jac), 0.90, 2) and same_at(max(jac), 0.97, 2),
        )
    ]
    details = {
        "available": True,
        "provenance_started_at": prov.get("started_at"),
        "recluster_version": (prov.get("environment") or {}).get("recluster"),
        "candidate_rule": rule,
        "n_runs": len(runs),
        "n_found": len(found),
        "jaccard": jac,
        "candidate_by": [r["match"]["candidate_by"] for r in runs],
    }
    return rows, details


# --------------------------------------------------------------------------- item 4


def _welch_se(a: np.ndarray, b: np.ndarray) -> float:
    return math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))


def item_pairing(kang: Kang, cpus: int) -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 6: paired vs unpaired composition tests and monocyte DE, labels as given."""
    counts, sample, condition = kang.counts, kang.sample, kang.condition
    thr = Thresholds()
    cmp_comp = Comparison(group_a="stim", group_b="ctrl", direction="different")
    pair = pair_key(
        sample,
        condition,
        counts.obs_names,
        "stim",
        "ctrl",
        SAMPLE_KEY,
        CONDITION,
        counts.obs,
    )
    assert pair is not None, "Kang should be a paired design"
    comp: dict[str, Any] = {}
    for cell_type in sorted(kang.labels.unique()):
        cells = counts.obs_names[(kang.labels == cell_type).to_numpy()]
        unp = evaluate_composition(cells, counts.obs_names, sample, condition, cmp_comp, thr)
        par = evaluate_composition(
            cells, counts.obs_names, sample, condition, cmp_comp, thr, pair=pair
        )
        table = pd.DataFrame(par["per_sample"])
        a = table[table["condition"] == "stim"].set_index("pair")
        b = table[table["condition"] == "ctrl"].set_index("pair")
        d = (a["logit"] - b["logit"].reindex(a.index)).to_numpy()
        se_unpaired = _welch_se(a["logit"].to_numpy(), b["logit"].to_numpy())
        se_paired = float(d.std(ddof=1) / math.sqrt(len(d)))
        lower = (a["proportion"] < b["proportion"].reindex(a.index)).to_numpy()
        comp[cell_type] = {
            "testable": bool(unp["evaluable"]),
            "design_paired": par["design"],
            "p_unpaired_welch": unp["p_logit_ttest"],
            "p_paired_t": par["p_logit_ttest"],
            "p_paired_signed_rank": par["p_mannwhitney"],
            "mean_proportion_stim": par["mean_proportion_a"],
            "mean_proportion_ctrl": par["mean_proportion_b"],
            "se_unpaired": se_unpaired,
            "se_paired": se_paired,
            "se_shrink_fold": se_unpaired / se_paired if se_paired > 0 else None,
            "n_donors_lower_in_stim": int(lower.sum()),
            "n_donors": int(len(lower)),
        }
    testable = {k: v for k, v in comp.items() if v["testable"]}
    shrink = [v["se_shrink_fold"] for v in testable.values()]
    p_unp = [v["p_unpaired_welch"] for v in testable.values()]
    paired_calls = [k for k, v in testable.items() if v["p_paired_t"] < thr.alpha]
    unpaired_calls = [k for k, v in testable.items() if v["p_unpaired_welch"] < thr.alpha]
    mono = comp[MONO_LABEL]

    # DE in the labelled monocytes over every gene with >= 10 counts, both designs
    say("  item 4: DESeq2 on labelled CD14+ monocytes, ~condition and ~pair + condition")
    cells = counts.obs_names[(kang.labels == MONO_LABEL).to_numpy()]
    pb, _ = pseudobulk_counts(counts, cells, sample, min_cells=thr.min_cells_per_sample)
    genes = [str(g) for g in pb.columns[pb.sum(axis=0) >= thr.pseudobulk_min_total_counts]]
    cmp_de = Comparison(group_a="stim", group_b="ctrl", direction="up_in_a")
    unp_df, unp_info = pseudobulk_de(
        counts, cells, sample, condition, cmp_de, genes, thr, n_cpus=cpus
    )
    par_df, par_info = pseudobulk_de(
        counts, cells, sample, condition, cmp_de, genes, thr, n_cpus=cpus, pair=pair
    )
    padj_u = unp_df["pb_padj_genomewide"]
    padj_p = par_df["pb_padj_genomewide"].reindex(padj_u.index)
    lfc_u = unp_df["pb_log2fc"]
    lfc_p = par_df["pb_log2fc"].reindex(lfc_u.index)
    called_u = padj_u < thr.alpha
    called_p = padj_p < thr.alpha
    same_sign = np.sign(lfc_u) == np.sign(lfc_p)
    kept = called_u & called_p & same_sign
    flips_called = called_u & called_p & ~same_sign
    flips_any = called_u & ~same_sign
    n_u = int(called_u.sum())
    share_kept = float(kept.sum() / n_u) if n_u else float("nan")
    de = {
        "n_genes_tested": len(genes),
        "design_unpaired": unp_info["design_formula"],
        "design_paired": par_info["design_formula"],
        "n_called_unpaired": n_u,
        "n_called_paired": int(called_p.sum()),
        "n_unpaired_calls_kept": int(kept.sum()),
        "share_unpaired_calls_kept": share_kept,
        "n_sign_flips_both_called": int(flips_called.sum()),
        "n_sign_differences_among_unpaired_calls": int(flips_any.sum()),
    }
    how_comp = (
        "evaluate_composition per cell_type (all cells as the run), sample = [replicate, label], "
        "without and with pair_key; SE of the logit difference: Welch vs per-donor differences"
    )
    rows = [
        row(
            "4 (RESULTS 6)",
            "Testable cell types; SE shrink, unpaired SE / paired SE",
            "7 types; 2- to 9-fold",
            f"{len(testable)} types; {min(shrink):.1f}- to {max(shrink):.1f}-fold",
            how_comp,
            len(testable) == 7 and round(min(shrink)) == 2 and round(max(shrink)) == 9,
        ),
        row(
            "4 (RESULTS 6)",
            "Unpaired Welch p over the testable types",
            "all >= 0.67, none called",
            f"min {min(p_unp):.2f}, {len(unpaired_calls)} called",
            "p_logit_ttest of the unpaired evaluate_composition",
            min(p_unp) >= 0.665 and not unpaired_calls,
        ),
        row(
            "4 (RESULTS 6)",
            "Paired calls; CD14+ Monocytes p, share stim vs ctrl, donors lower in stim",
            "CD14+ Monocytes only; p = 0.035, 23.8 % vs 25.4 %, 7 of 8",
            f"{', '.join(paired_calls) or 'none'}; p = {mono['p_paired_t']:.3f}, "
            f"{100 * mono['mean_proportion_stim']:.1f} % vs "
            f"{100 * mono['mean_proportion_ctrl']:.1f} %, "
            f"{mono['n_donors_lower_in_stim']} of {mono['n_donors']}",
            "paired evaluate_composition (paired t on per-donor logit differences)",
            paired_calls == [MONO_LABEL]
            and same_at(mono["p_paired_t"], 0.035, 3)
            and same_at(100 * mono["mean_proportion_stim"], 23.8, 1)
            and same_at(100 * mono["mean_proportion_ctrl"], 25.4, 1)
            and mono["n_donors_lower_in_stim"] == 7,
        ),
        row(
            "4 (RESULTS 6)",
            "Monocyte DE: share of unpaired calls kept under ~pair + condition; sign flips",
            "95 %; 0",
            f"{100 * share_kept:.1f} % ({int(kept.sum())} of {n_u}); "
            f"{int(flips_called.sum())} flips",
            f"pseudobulk_de on cells labelled CD14+ Monocytes, {len(genes)} genes with >= 10 "
            "counts; call = DESeq2 padj < 0.05; kept = also called with the same sign",
            same_at(100 * share_kept, 95, 0) and int(flips_called.sum()) == 0,
        ),
    ]
    return rows, {"composition": comp, "de_monocytes": de}


# --------------------------------------------------------------------------- item 5


def item_doublets(data_dir: Path) -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 8: expected doublet rate vs demuxlet, and the planted PBMC3k doublets."""
    grid = GridConfig()
    say("  item 5: Scrublet on the two GEO lanes")
    geo = load_batch2(data_dir / "geo")
    pipe = Pipeline(geo, grid, sample=geo.obs["lane"], mito_prefix=None)
    stage = pipe.qc_stage(grid.baseline)
    obs = stage.adata.obs
    lanes: dict[str, Any] = {}
    for lane in sorted(obs[SAMPLE_COL].unique()):
        sel = obs[obs[SAMPLE_COL] == lane]
        n = len(sel)
        expected = min(0.30, grid.scrublet_expected_rate_per_1000 * n / 1000.0)
        agree = _agreement(sel)
        lanes[lane] = {
            "n_scored": n,
            "expected_rate": expected,
            "demuxlet_doublet_rate": agree["demuxlet_doublet_rate"],
            "difference_pp": 100 * (expected - agree["demuxlet_doublet_rate"]),
            "called_rate_all_scored": float(sel["predicted_doublet"].mean()),
            "agreement": agree,
        }
    worst = max(abs(v["difference_pp"]) for v in lanes.values())

    say("  item 5: Scrublet on the planted PBMC3k doublets")
    planted = load_adata(data_dir / "planted" / "pbmc3k_doublets.h5ad")
    counts, _ = extract_counts(planted)
    ppipe = Pipeline(counts, grid, mito_prefix=detect_mito_prefix(counts.var_names))
    pstage = ppipe.qc_stage(grid.baseline)
    pobs = pstage.adata.obs
    injected = pobs["injected_doublet"].to_numpy(dtype=bool)
    score = pobs["doublet_score"].to_numpy(dtype=float)
    k = int(injected.sum())
    top = np.argsort(-score, kind="stable")[:k]
    oracle_recall = float(injected[top].mean())
    median_injected = float(np.median(score[injected]))
    pdet = {
        "n_cells_after_qc": int(len(pobs)),
        "n_injected_after_qc": k,
        "median_score_injected": median_injected,
        "median_score_other": float(np.median(score[~injected])),
        "oracle_top_k_recall": oracle_recall,
        "pipeline_recall": float(pobs["predicted_doublet"].to_numpy(dtype=bool)[injected].mean()),
        "auroc_of_score": auroc_two_sets(score, injected),
    }
    lane_text = "; ".join(
        f"{lane} {100 * v['expected_rate']:.1f} % vs {100 * v['demuxlet_doublet_rate']:.1f} %"
        for lane, v in lanes.items()
    )
    rows = [
        row(
            "5 (RESULTS 8)",
            "Expected doublet rate vs demuxlet, per GEO lane (ambiguous excluded)",
            "within 0.7 percentage points",
            f"{lane_text}; largest gap {worst:.2f} pp",
            "load_batch2 + Pipeline.qc_stage at baseline QC, Scrublet per lane; expected = 0.8 % "
            "per 1,000 scored cells of the lane; demuxlet rate over singlets + doublets",
            worst <= 0.75,
        ),
        row(
            "5 (RESULTS 8)",
            "Planted PBMC3k doublets: median Scrublet score of injected cells",
            "0.10",
            f"{median_injected:.3f}",
            "Pipeline.qc_stage at baseline QC on data/planted/pbmc3k_doublets.h5ad (one subset)",
            same_at(median_injected, 0.10, 2),
        ),
        row(
            "5 (RESULTS 8)",
            "Share of injected cells caught by an oracle cutoff at the true count",
            "28 %",
            f"{100 * oracle_recall:.1f} % (top {k} scores)",
            "injected share among the k highest scores, k = injected cells passing QC",
            same_at(100 * oracle_recall, 28, 0),
        ),
    ]
    return rows, {"geo_lanes": lanes, "planted_pbmc3k": pdet}


# --------------------------------------------------------------------------- item 6


def _log2_cpm_by(counts: AnnData, cells: pd.Index, group: pd.Series, gene: str) -> pd.Series:
    pb, _ = pseudobulk_counts(counts, cells, group, min_cells=1)
    return np.log2(pb[gene] / pb.sum(axis=1) * 1e6 + 1.0)


def item_time_gene(kang: Kang, cpus: int) -> tuple[list[Row], dict[str, Any]]:
    """RESULTS item 12: RPS26 vs SQSTM1, and the planted SQSTM1 effect at half-amplitude 0.4."""
    counts = kang.counts
    ctrl = (kang.condition == "ctrl").to_numpy()
    ctrl_cells = counts.obs_names[ctrl]
    ctrl_mono = counts.obs_names[ctrl & (kang.labels == MONO_LABEL).to_numpy()]
    spread: dict[str, dict[str, float]] = {}
    for gene in ("RPS26", "SQSTM1"):
        spread[gene] = {}
        for name, cells in (("control cells", ctrl_cells), ("control monocytes", ctrl_mono)):
            v = _log2_cpm_by(counts, cells, kang.donor, gene)
            spread[gene][name] = float(v.max() - v.min())

    thr = Thresholds()
    cmp = Comparison(
        group_a=kang_module.TIME_CLAIM_GROUP_A,
        group_b=kang_module.TIME_CLAIM_GROUP_B,
        direction="up_in_a",
    )
    gene = kang_module.TIME_CLAIM_GENE
    planted: dict[str, Any] = {}
    original = kang_module.TIME_SIGNAL
    for amp in (0.4, 0.6):
        say(f"  item 6: planting {gene} at half-amplitude {amp}")
        signal = dict(original)
        signal[gene] = (original[gene][0], amp)
        kang_module.TIME_SIGNAL = signal
        try:
            p_adata = kang_module.plant_time_signal(kang.adata, seed=0)
        finally:
            kang_module.TIME_SIGNAL = original
        p_counts, _ = extract_counts(p_adata)
        obs = p_counts.obs
        cells = p_counts.obs_names[(obs["cell_type"].astype(str) == MONO_LABEL).to_numpy()]
        df, info = pseudobulk_de(
            p_counts,
            cells,
            obs["replicate"].astype(str),
            obs["fake_group"].astype(str),
            cmp,
            [gene],
            thr,
            n_cpus=cpus,
        )
        planted[str(amp)] = {
            "design": info["design_formula"],
            "n_samples_a": info["n_samples_a"],
            "n_samples_b": info["n_samples_b"],
            "log2fc": float(df.loc[gene, "pb_log2fc"]),
            "p": float(df.loc[gene, "pb_p"]),
            "padj_claimed_genes": float(df.loc[gene, "pb_p"]),  # BH over one gene is the p
            "padj_genomewide": float(df.loc[gene, "pb_padj_genomewide"]),
        }
    low = planted["0.4"]
    rows = [
        row(
            "6 (RESULTS 12)",
            "RPS26 donor-to-donor log2 CPM range, control cells (control monocytes)",
            "about 2.6 (kang.py comment: 1.7 in monocytes)",
            f"{spread['RPS26']['control cells']:.2f} ({spread['RPS26']['control monocytes']:.2f})",
            "log2(CPM + 1) of per-donor pseudobulks of the ctrl cells, max - min over 8 donors",
            same_at(spread["RPS26"]["control cells"], 2.6, 1),
        ),
        row(
            "6 (RESULTS 12)",
            "SQSTM1 donor-to-donor log2 CPM range, control monocytes (control cells)",
            "about 0.4 (kang.py comment: 0.6 over all control cells)",
            f"{spread['SQSTM1']['control monocytes']:.2f} "
            f"({spread['SQSTM1']['control cells']:.2f})",
            "as above",
            same_at(spread["SQSTM1"]["control monocytes"], 0.4, 1),
        ),
        row(
            "6 (RESULTS 12)",
            "Planted SQSTM1 at half-amplitude 0.4: unadjusted pseudobulk padj",
            "0.065",
            f"{low['padj_claimed_genes']:.3f} (log2FC {low['log2fc']:.2f}; genome-wide padj "
            f"{low['padj_genomewide']:.3f}); at 0.6: {planted['0.6']['padj_claimed_genes']:.4f}"
            f" (log2FC {planted['0.6']['log2fc']:.2f})",
            "plant_time_signal(seed=0) with TIME_SIGNAL['SQSTM1'] half-amplitude replaced; "
            "pseudobulk_de on cells labelled CD14+ Monocytes, sample = replicate, B vs A (4 vs 4), "
            "~condition; padj over the one claimed gene",
            same_at(low["padj_claimed_genes"], 0.065, 3),
        ),
    ]
    return rows, {"donor_log2cpm_range": spread, "planted_sqstm1": planted}


# --------------------------------------------------------------------------- item 7


def item_clock_panel(kang: Kang) -> tuple[list[Row], dict[str, Any]]:
    """HOW_IT_WORKS section 8: Mode B on the real IFN-beta comparison."""
    counts = kang.counts
    thr = Thresholds()
    unit, unit_col, _ = exchangeable_units(
        counts.obs,
        counts.obs_names,
        kang.sample,
        kang.condition,
        "stim",
        "ctrl",
        SAMPLE_KEY,
        CONDITION,
    )
    group_by_unit = kang.condition.loc[unit.index].groupby(unit.to_numpy()).first()
    panel = panel_check(counts, unit, group_by_unit, "stim", "ctrl", thr, unit_column=unit_col)
    padj = {r["gene"]: float(r["padj"]) for r in panel["per_gene"]}
    cited = {"DBP": 0.012, "CRY1": 0.037, "PER1": 0.051}
    how = (
        "exchangeable_units + panel_check on all Kang cells, sample = [replicate, label], "
        f"stim vs ctrl ({panel['n_samples_a']} vs {panel['n_samples_b']} units)"
    )
    rows = [
        row(
            "7 (HOW_IT_WORKS 8)",
            f"Clock-gene per-gene BH padj, {g}",
            f"{c:.3f}",
            f"{padj[g]:.3f}",
            how,
            same_at(padj[g], c, 3),
        )
        for g, c in cited.items()
    ]
    rows.append(
        row(
            "7 (HOW_IT_WORKS 8)",
            "Clock-state index permutation p",
            "0.38",
            f"{panel['index_p_permutation']:.3f} ({panel['index_n_permutations']} "
            f"{'exact' if panel['index_exact'] else 'sampled'} permutations)",
            how,
            same_at(float(panel["index_p_permutation"]), 0.38, 2),
        )
    )
    return rows, {"unit_column": unit_col, "panel": panel}


# --------------------------------------------------------------------------- output


def write_markdown(rows: list[Row], meta: dict[str, Any], path: Path) -> None:
    lines = [
        "# Design measurements",
        "",
        'Numbers cited in `bench/RESULTS.md` ("What changed between the runs") and '
        "`docs/HOW_IT_WORKS.md` (section 8), recomputed from `data/` by "
        "`scripts/design_measurements.py` with the package's own functions. Unrounded values "
        'and per-item details are in `design_measurements.json`. "Agrees" means the measured '
        "value rounds to the cited one at the cited precision (or lies inside a cited bound).",
        "",
        f"Recluster {meta['recluster_version']}, run {meta['finished_at']}, "
        f"{meta['runtime_s'] / 60:.1f} min. The cluster-level doublet mixture analysis "
        "(RESULTS item 9) is not part of this script.",
        "",
        "| Item | Quantity | Cited | Measured | Agrees | How computed |",
        "|---|---|---|---|---|---|",
    ]
    for r in rows:
        agrees = "–" if r["agrees"] is None else ("yes" if r["agrees"] else "**no**")
        cells = [r["item"], r["quantity"], r["cited"], r["measured"], agrees, r["how"]]
        lines.append("| " + " | ".join(c.replace("|", "\\|") for c in cells) + " |")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--bench", type=Path, default=ROOT / "bench")
    ap.add_argument("--out", type=Path, default=None, help="output folder (default: --bench)")
    ap.add_argument("--only", default="1,2,3,4,5,6,7", help="comma-separated item numbers")
    ap.add_argument("--cpus", type=int, default=3, help="threads for pydeseq2")
    args = ap.parse_args()
    out_dir = args.out or args.bench
    wanted = {int(x) for x in args.only.split(",") if x.strip()}

    t0 = time.perf_counter()
    kang: Kang | None = None

    def need_kang() -> Kang:
        nonlocal kang
        if kang is None:
            say("loading Kang 2018")
            kang = Kang(args.data)
        return kang

    items: dict[int, Callable[[], tuple[list[Row], dict[str, Any]]]] = {
        1: lambda: item_qc_cap(need_kang(), args.cpus),
        2: item_simulation,
        3: lambda: item_cd4t_jaccard(args.bench),
        4: lambda: item_pairing(need_kang(), args.cpus),
        5: lambda: item_doublets(args.data),
        6: lambda: item_time_gene(need_kang(), args.cpus),
        7: lambda: item_clock_panel(need_kang()),
    }
    rows: list[Row] = []
    details: dict[str, Any] = {}
    timings: dict[str, float] = {}
    for i, fn in items.items():
        if i not in wanted:
            continue
        say(f"item {i}")
        t = time.perf_counter()
        r, d = fn()
        timings[str(i)] = round(time.perf_counter() - t, 1)
        rows += r
        details[str(i)] = d
        for x in r:
            say(f"  {x['quantity']}: cited {x['cited']}; measured {x['measured']}")

    meta = {
        "recluster_version": recluster.__version__,
        "finished_at": time.strftime("%Y-%m-%d %H:%M"),
        "runtime_s": round(time.perf_counter() - t0, 1),
        "item_runtime_s": timings,
        "items_run": sorted(wanted),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"meta": meta, "rows": rows, "details": details}
    (out_dir / "design_measurements.json").write_text(
        json.dumps(payload, indent=2, default=_json_default) + "\n"
    )
    write_markdown(rows, meta, out_dir / "design_measurements.md")
    say(f"wrote {out_dir / 'design_measurements.json'} and design_measurements.md")


def _json_default(x: Any) -> Any:
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, (pd.Index, np.ndarray)):
        return list(x)
    return str(x)


if __name__ == "__main__":
    main()
