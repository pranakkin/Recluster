"""Orchestration: load -> grid -> run -> match -> evaluate -> verdict -> report."""

from __future__ import annotations

import csv
import logging
import re
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd

from recluster.config import GridConfig, Thresholds
from recluster.data import DataError, detect_mito_prefix, extract_counts, load_adata, sha256_file
from recluster.evaluators.composition import evaluate_composition
from recluster.evaluators.de import evaluate_de
from recluster.evaluators.population import evaluate_population
from recluster.grid import build_grid, describe_grid
from recluster.matching import Match, lognorm_matrix, marker_score, match_population
from recluster.pipeline import Pipeline, RunResult
from recluster.provenance import dump_json, environment, round_floats, timestamp
from recluster.report.html import write_report
from recluster.report.methods import methods_paragraph
from recluster.results import ClaimResult, RunEvaluation
from recluster.spec import ClaimSpec, make_sample_key
from recluster.verdict import compute_verdict

log = logging.getLogger(__name__)
Progress = Callable[[str], None]


def _stderr(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def slugify(name: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "-", name).strip("-").lower()
    return s[:60] or "claim"


def run_claim(
    data_path: str | Path,
    spec: ClaimSpec,
    *,
    out_dir: str | Path = "results",
    grid: GridConfig | None = None,
    thresholds: Thresholds | None = None,
    audit: bool = False,
    progress: Progress | None = _stderr,
    write: bool = True,
) -> ClaimResult:
    """Stress-test one claim and (optionally) write report.html, provenance.json, runs.csv."""
    t_start = time.perf_counter()
    started_at = timestamp()
    say = progress or (lambda _m: None)
    grid = grid or GridConfig()
    thr = thresholds or Thresholds()
    data_path = Path(data_path)

    # ---------------------------------------------------------------- data
    say(f"loading {data_path}")
    adata = load_adata(data_path)
    spec.check_against(adata, audit=audit)
    counts, counts_source = extract_counts(adata)
    if audit and counts_source == ".X":
        raise DataError(
            "audit mode expects a processed object whose raw counts live in layers['counts'] "
            "or .raw; "
            ".X itself looks like raw counts, so use `recluster run` instead"
        )
    sha = sha256_file(data_path)
    cols = spec.columns
    obs = counts.obs
    sample = make_sample_key(obs, cols.sample_key) if cols.sample_key else None
    condition = obs[cols.condition].astype(str) if cols.condition else None
    batch = obs[cols.batch].astype(str) if cols.batch else None
    labels = obs[cols.labels].astype(str) if cols.labels else None
    mito_prefix = detect_mito_prefix(counts.var_names, grid.mito_prefix)

    data_summary: dict[str, Any] = {
        "path": str(data_path),
        "path_name": data_path.name,
        "sha256": sha,
        "n_cells": int(counts.n_obs),
        "n_genes": int(counts.n_vars),
        "counts_source": counts_source,
        "mito_prefix": mito_prefix,
        "n_samples": int(sample.nunique()) if sample is not None else None,
        "conditions": sorted(condition.unique().tolist()) if condition is not None else None,
        "n_batches": int(batch.nunique()) if batch is not None else None,
    }

    # ------------------------------------------------------- expression once
    tgt = spec.target
    needed = list(dict.fromkeys(tgt.markers + tgt.negative_markers + spec.genes))
    lognorm = lognorm_matrix(counts, needed)
    score = marker_score(lognorm, tgt.markers, tgt.negative_markers) if tgt.markers else None

    reference: pd.Index | None = None
    reference_source = "none"
    if tgt.label_value is not None and labels is not None:
        reference = counts.obs_names[(labels == tgt.label_value).to_numpy()]
        reference_source = (
            f"cells labelled `{tgt.label_value}` in `{cols.labels}` ({len(reference)} cells)"
        )

    # ---------------------------------------------------------------- grid
    runs_params = build_grid(grid, has_batch=batch is not None, has_mito=mito_prefix is not None)
    grid_counts = describe_grid(runs_params)
    pipe = Pipeline(counts, grid, sample=sample, batch=batch, mito_prefix=mito_prefix)
    say(
        f"grid: {grid_counts['total']} runs ({grid_counts['one_at_a_time']} one-at-a-time, "
        f"{grid_counts['random']} random)"
    )

    evaluations: list[RunEvaluation] = []
    base_params = runs_params[0]
    for i, params in enumerate(runs_params):
        t0 = time.perf_counter()
        rr = pipe.run(params)
        m = match_population(
            rr.clusters,
            score=score,
            reference=reference,
            sample=sample,
            min_cells=tgt.min_cells,
            thr=thr,
        )
        if (
            reference is None
            and m.found
            and m.marker_auroc is not None
            and m.marker_auroc >= thr.marker_auroc_min
        ):
            reference = m.cells
            where = "the baseline run" if i == 0 else f"run {i}"
            reference_source = f"population found in {where} ({len(reference)} cells)"
            m = match_population(
                rr.clusters,
                score=score,
                reference=reference,
                sample=sample,
                min_cells=tgt.min_cells,
                thr=thr,
            )
        ev = _evaluate(spec, m, rr, counts, lognorm, sample, condition, thr, grid)
        match_info = m.to_dict()
        match_info.update(_doublet_summary(rr.doublet_scores, m.cells))
        label = params.label(base_params)
        evaluations.append(
            RunEvaluation(
                index=i,
                label=label,
                params=params.model_dump(),
                is_baseline=(i == 0),
                is_one_at_a_time=(len(params.diff(base_params)) == 1),
                n_cells=rr.n_cells,
                n_removed_qc=rr.n_removed_qc,
                n_removed_doublet=rr.n_removed_doublet,
                n_clusters=rr.n_clusters,
                runtime_s=round(time.perf_counter() - t0, 3),
                warnings=rr.warnings + list(ev.get("warnings", [])),
                match=match_info,
                holds=bool(ev["holds"]),
                holds_cell_level=ev.get("holds_cell_level"),
                replicate_ok=ev.get("replicate_ok"),
                recovered=ev.get("recovered"),
                key_stats=ev.get("key_stats", {}),
                details={
                    k: v
                    for k, v in ev.items()
                    if k
                    not in (
                        "holds",
                        "holds_cell_level",
                        "replicate_ok",
                        "recovered",
                        "key_stats",
                        "warnings",
                    )
                },
            )
        )
        say(
            f"run {i:>2}/{len(runs_params) - 1}  {label:<40.40s} cells={rr.n_cells:<6} "
            f"k={rr.n_clusters:<3} pop={m.n_cells:<5} "
            f"{'HELD' if ev['holds'] else 'failed':<6} {time.perf_counter() - t0:5.1f}s"
        )

    # ------------------------------------------------------------ as labelled
    as_labelled = None
    if tgt.label_value is not None and labels is not None:
        as_labelled = _as_labelled(spec, counts, lognorm, sample, condition, thr, grid, labels)

    verdict = compute_verdict(evaluations, base_params, thr, spec.claim_type)
    methods = methods_paragraph(
        spec, grid, thr, data_summary, len(runs_params), grid_counts, reference_source
    )
    total = time.perf_counter() - t_start
    provenance = {
        "started_at": started_at,
        "finished_at": timestamp(),
        "total_runtime_s": round(total, 3),
        "mode": "audit" if audit else "run",
        "environment": environment(),
        "input": data_summary,
        "claim": spec.model_dump(exclude_none=True),
        "grid": grid.model_dump(),
        "thresholds": thr.model_dump(),
        "grid_counts": grid_counts,
        "reference_source": reference_source,
        "verdict": verdict.to_dict(),
        "runs": [round_floats(e.to_dict()) for e in evaluations],
        "as_labelled": round_floats(as_labelled) if as_labelled else None,
    }
    result = ClaimResult(
        spec=spec,
        mode="audit" if audit else "run",
        verdict=verdict,
        runs=evaluations,
        as_labelled=as_labelled,
        reference_source=reference_source,
        data_summary=data_summary,
        grid=grid,
        thresholds=thr,
        provenance=provenance,
        methods=methods,
        total_runtime_s=total,
    )
    if write:
        out = Path(out_dir) / slugify(spec.name)
        out.mkdir(parents=True, exist_ok=True)
        write_report(result, out / "report.html")
        dump_json(provenance, out / "provenance.json")
        (out / "claim.yaml").write_text(spec.to_yaml())
        _write_runs_csv(evaluations, out / "runs.csv")
        say(f"{verdict.label}: {verdict.summary}")
        say(f"wrote {out / 'report.html'}")
    return result


def _evaluate(
    spec: ClaimSpec,
    m: Match,
    rr: RunResult,
    counts: Any,
    lognorm: pd.DataFrame,
    sample: pd.Series | None,
    condition: pd.Series | None,
    thr: Thresholds,
    grid: GridConfig,
) -> dict[str, Any]:
    if spec.claim_type == "population_exists":
        return evaluate_population(m, thr)
    assert spec.comparison is not None and sample is not None and condition is not None
    if spec.claim_type == "composition_shift":
        return evaluate_composition(
            m.cells, rr.clusters.index, sample, condition, spec.comparison, thr
        )
    if not m.found:
        return {
            "holds": False,
            "holds_cell_level": False,
            "replicate_ok": False,
            "reasons": [f"population not found ({m.n_cells} cells)"],
            "key_stats": {"pb_log2fc": None, "pb_padj": None, "cell_padj": None},
            "per_gene": [],
        }
    return evaluate_de(
        cells=m.cells,
        counts=counts,
        lognorm=lognorm,
        sample=sample,
        condition=condition,
        comparison=spec.comparison,
        genes=spec.genes,
        require=spec.require,
        thr=thr,
        n_cpus=grid.n_cpus,
    )


def _as_labelled(
    spec: ClaimSpec,
    counts: Any,
    lognorm: pd.DataFrame,
    sample: pd.Series | None,
    condition: pd.Series | None,
    thr: Thresholds,
    grid: GridConfig,
    labels: pd.Series,
) -> dict[str, Any]:
    cells = counts.obs_names[(labels == spec.target.label_value).to_numpy()]
    all_cells = counts.obs_names
    if spec.claim_type == "population_exists":
        out: dict[str, Any] = {"n_cells": int(len(cells))}
        if sample is not None and len(cells):
            per = sample.loc[cells].value_counts()
            out["n_samples_present"] = int((per >= thr.min_cells_per_sample).sum())
            out["max_sample_share"] = float(per.max() / len(cells))
        return out
    assert spec.comparison is not None and sample is not None and condition is not None
    if spec.claim_type == "composition_shift":
        return evaluate_composition(cells, all_cells, sample, condition, spec.comparison, thr)
    return evaluate_de(
        cells=cells,
        counts=counts,
        lognorm=lognorm,
        sample=sample,
        condition=condition,
        comparison=spec.comparison,
        genes=spec.genes,
        require=spec.require,
        thr=thr,
        n_cpus=grid.n_cpus,
    )


def _doublet_summary(scores: pd.Series, cells: pd.Index) -> dict[str, float | None]:
    """Median Scrublet score of the matched population vs all kept cells (informational)."""
    if scores.isna().all():
        return {"doublet_score_median": None, "doublet_score_median_all": None}
    pop = scores.loc[scores.index.intersection(cells)]
    return {
        "doublet_score_median": None if pop.empty else float(pop.median()),
        "doublet_score_median_all": float(scores.median()),
    }


def _write_runs_csv(evaluations: list[RunEvaluation], path: Path) -> None:
    if not evaluations:
        return
    param_keys = list(evaluations[0].params.keys())
    stat_keys = list(evaluations[0].key_stats.keys())
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "index",
                "label",
                *param_keys,
                "n_cells",
                "n_clusters",
                "population_cells",
                "jaccard",
                "marker_auroc",
                *stat_keys,
                "holds",
                "holds_cell_level",
                "runtime_s",
            ]
        )
        for e in evaluations:
            w.writerow(
                [
                    e.index,
                    e.label,
                    *[e.params[k] for k in param_keys],
                    e.n_cells,
                    e.n_clusters,
                    e.match.get("n_cells"),
                    e.match.get("jaccard"),
                    e.match.get("marker_auroc"),
                    *[e.key_stats.get(k) for k in stat_keys],
                    e.holds,
                    e.holds_cell_level,
                    e.runtime_s,
                ]
            )
