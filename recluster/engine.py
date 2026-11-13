"""Orchestration: load -> grid -> run -> match -> evaluate -> verdict -> report."""

from __future__ import annotations

import csv
import multiprocessing
import os
import re
import sys
import time
from collections.abc import Callable, Iterator
from concurrent.futures import Future, ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from recluster.config import GridConfig, PipelineParams, Thresholds
from recluster.data import (
    AUDIT_COUNTS_IN_X,
    DataError,
    detect_mito_prefix,
    extract_counts,
    load_adata,
    sha256_file,
)
from recluster.evaluators.composition import evaluate_composition
from recluster.evaluators.de import evaluate_de
from recluster.evaluators.doublet_cluster import doublet_cluster_check
from recluster.evaluators.pairing import pair_key
from recluster.evaluators.population import evaluate_population
from recluster.grid import build_grid, describe_grid
from recluster.matching import Match, lognorm_matrix, marker_score, match_population
from recluster.pipeline import Pipeline, RunResult
from recluster.provenance import dump_json, environment, round_floats, timestamp
from recluster.report.html import write_report
from recluster.report.methods import methods_paragraph
from recluster.results import ClaimResult, RunEvaluation
from recluster.spec import ClaimSpec, make_sample_key
from recluster.time_confound import assess_claim, load_rhythmic_table
from recluster.verdict import apply_time_check, compute_verdict

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
    jobs: int = 1,
) -> ClaimResult:
    """Stress-test one claim and (optionally) write report.html, provenance.json, runs.csv.

    With ``jobs`` > 1 the grid points run in that many worker processes. Results are the same
    as with ``jobs=1``: every run is seeded and independent of the others, matching stays in
    grid order in this process, and only timings differ.
    """
    if jobs < 1:
        raise ValueError(f"jobs must be at least 1, got {jobs}")
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
        raise DataError(AUDIT_COUNTS_IN_X)
    sha = sha256_file(data_path)
    cols = spec.columns
    obs = counts.obs
    sample = make_sample_key(obs, cols.sample_key) if cols.sample_key else None
    condition = obs[cols.condition].astype(str) if cols.condition else None
    batch = obs[cols.batch].astype(str) if cols.batch else None
    labels = obs[cols.labels].astype(str) if cols.labels else None
    mito_prefix = detect_mito_prefix(counts.var_names, grid.mito_prefix)

    pair: pd.Series | None = None
    if spec.comparison is not None and sample is not None and condition is not None:
        pair = pair_key(
            sample,
            condition,
            counts.obs_names,
            spec.comparison.group_a,
            spec.comparison.group_b,
            cols.sample_key,
            cols.condition,
            obs,
        )

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
        "paired_design": pair is not None,
        "n_pairs": int(pair.nunique()) if pair is not None else None,
    }

    # ------------------------------------------------------- expression once
    tgt = spec.target
    needed = list(dict.fromkeys(tgt.markers + tgt.negative_markers + spec.genes))
    lognorm = lognorm_matrix(counts, needed)
    score = (
        marker_score(lognorm, tgt.markers, tgt.negative_markers, logic=tgt.marker_logic)
        if tgt.markers
        else None
    )
    marker_frame = lognorm[tgt.markers] if tgt.markers else None

    reference: pd.Index | None = None
    reference_source = "none"
    if tgt.label_value is not None and labels is not None:
        reference = counts.obs_names[(labels == tgt.label_value).to_numpy()]
        reference_source = (
            f"cells labelled `{tgt.label_value}` in `{cols.labels}` ({len(reference)} cells)"
        )

    # ---------------------------------------------------------------- grid
    pipe = Pipeline(counts, grid, sample=sample, batch=batch, mito_prefix=mito_prefix)
    runs_params, merged_runs = _dedupe_after_capping(
        build_grid(grid, has_batch=batch is not None, has_mito=mito_prefix is not None), pipe
    )
    grid_counts = describe_grid(runs_params)
    grid_counts["merged_after_capping"] = len(merged_runs)
    say(
        f"grid: {grid_counts['total']} runs ({grid_counts['one_at_a_time']} one-at-a-time, "
        f"{grid_counts['random']} random"
        + (f", {len(merged_runs)} dropped as duplicates after QC capping" if merged_runs else "")
        + ")"
        + ("; paired design" if pair is not None else "")
    )
    reference_is_labelled = reference is not None

    def match(rr: RunResult, ref: pd.Index | None) -> Match:
        return match_population(
            rr.clusters,
            score=score,
            reference=ref,
            sample=sample,
            min_cells=tgt.min_cells,
            thr=thr,
            marker_frame=marker_frame,
            logic=tgt.marker_logic,
            reference_is_labelled=reference_is_labelled,
        )

    def match_run(i: int, rr: RunResult) -> Match:
        """Match run ``i``; must be called in grid order, since the first qualifying marker
        match becomes the reference for every later run."""
        nonlocal reference, reference_source
        m = match(rr, reference)
        if (
            reference is None
            and m.found
            and m.marker_auroc is not None
            and m.marker_auroc >= thr.marker_auroc_min
        ):
            reference = m.cells
            where = "the baseline run" if i == 0 else f"run {i}"
            reference_source = f"population found in {where} ({len(reference)} cells)"
            m = match(rr, reference)
        return m

    evaluations: list[RunEvaluation] = []
    base_params = runs_params[0]
    baseline_cells: pd.Index | None = None
    baseline_run_cells: pd.Index = counts.obs_names

    def finish(
        i: int,
        rr: RunResult,
        m: Match,
        ev: dict[str, Any],
        doublet_info: dict[str, Any] | None,
        runtime_s: float,
    ) -> None:
        nonlocal baseline_cells, baseline_run_cells
        if i == 0:
            baseline_cells = m.cells if m.found else None
            baseline_run_cells = rr.clusters.index
        e = _run_evaluation(i, runs_params[i], base_params, rr, m, ev, doublet_info, runtime_s)
        evaluations.append(e)
        say(
            f"run {i:>2}/{len(runs_params) - 1}  {e.label:<40.40s} cells={rr.n_cells:<6} "
            f"k={rr.n_clusters:<3} pop={m.n_cells:<5} "
            f"{'HELD' if ev['holds'] else 'failed':<6} {runtime_s:5.1f}s"
        )

    workers = min(jobs, len(runs_params))
    if workers == 1:
        inputs = _EvalInputs(spec, counts, lognorm, sample, condition, thr, grid.n_cpus, pair)
        with _single_threaded_blas():
            for i, params in enumerate(runs_params):
                t0 = time.perf_counter()
                rr = pipe.run(params)
                m = match_run(i, rr)
                ev, doublet_info = _evaluate_run(inputs, rr.clusters, m)
                finish(i, rr, m, ev, doublet_info, time.perf_counter() - t0)
    else:
        say(f"running on {workers} worker processes")
        initargs = (counts, grid, sample, batch, mito_prefix, spec, lognorm, condition, thr, pair)
        with _worker_pool(workers, initargs) as pool:
            run_futures = [pool.submit(_pipeline_task, p) for p in runs_params]
            # each run is matched as soon as it and all earlier runs are done; its evaluation
            # then goes back to the pool while later runs are still clustering
            matched: list[tuple[RunResult, Match, float, Future[_EvalOutput]]] = []
            for i, fut in enumerate(run_futures):
                rr = fut.result()
                t0 = time.perf_counter()
                m = match_run(i, rr)
                match_s = time.perf_counter() - t0
                matched.append((rr, m, match_s, pool.submit(_evaluate_task, rr.clusters, m)))
            for i, (rr, m, match_s, ev_future) in enumerate(matched):
                ev, doublet_info, eval_s = ev_future.result()
                finish(i, rr, m, ev, doublet_info, rr.runtime_s + match_s + eval_s)

    # ------------------------------------------------------------ as labelled
    as_labelled = None
    if tgt.label_value is not None and labels is not None:
        as_labelled = _as_labelled(
            spec, counts, lognorm, sample, condition, thr, grid, labels, pair
        )

    # ------------------------------------------------------ time of day
    time_check = None
    if evaluations and spec.claim_type != "population_exists":
        data_summary["baseline_design"] = evaluations[0].details.get("design")
    if spec.comparison is not None and sample is not None and condition is not None:
        say("time-of-day check")
        tc_cells = baseline_cells
        tc_run_cells = baseline_run_cells
        held_unadjusted: bool | None = evaluations[0].holds if baseline_cells is not None else None
        if tc_cells is None and tgt.label_value is not None and labels is not None:
            tc_cells = counts.obs_names[(labels == tgt.label_value).to_numpy()]
            tc_run_cells = counts.obs_names
        time_check = assess_claim(
            spec,
            counts,
            sample,
            condition,
            obs,
            tc_cells,
            tc_run_cells,
            thr,
            n_cpus=grid.n_cpus,
            rhythmic_table=load_rhythmic_table(),
            pair=pair,
            held_unadjusted=held_unadjusted,
        ).to_dict()

    verdict = apply_time_check(
        compute_verdict(evaluations, base_params, thr, spec.claim_type), time_check
    )
    methods = methods_paragraph(
        spec, grid, thr, data_summary, len(runs_params), grid_counts, reference_source
    )
    total = time.perf_counter() - t_start
    provenance = {
        "started_at": started_at,
        "finished_at": timestamp(),
        "total_runtime_s": round(total, 3),
        "mode": "audit" if audit else "run",
        "execution": {"jobs": jobs, "workers": workers},
        "environment": environment(),
        "input": data_summary,
        "claim": spec.model_dump(exclude_none=True),
        "grid": grid.model_dump(),
        "thresholds": thr.model_dump(),
        "grid_counts": grid_counts,
        "merged_runs": merged_runs,
        "reference_source": reference_source,
        "verdict": verdict.to_dict(),
        "time_check": round_floats(time_check) if time_check else None,
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
        time_check=time_check,
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


def _dedupe_after_capping(
    runs_params: list[PipelineParams], pipe: Pipeline
) -> tuple[list[PipelineParams], list[dict[str, Any]]]:
    """Drop grid points that become identical once the QC cell-loss cap is applied.

    Two requested thresholds that both cap to the same effective value would give bitwise
    identical runs and count one configuration twice in the survival rate. The first point is
    kept (the baseline always is) and each dropped point is recorded for the provenance.
    """
    base = runs_params[0]
    kept: list[PipelineParams] = []
    seen: dict[PipelineParams, int] = {}
    merged: list[dict[str, Any]] = []
    for p in runs_params:
        eff = pipe.effective_params(p)
        if eff in seen:
            merged.append(
                {
                    "dropped": p.label(base),
                    "same_as": kept[seen[eff]].label(base),
                    "effective": eff.diff(base) or {"note": "same as baseline"},
                }
            )
            continue
        seen[eff] = len(kept)
        kept.append(p)
    return kept, merged


@dataclass
class _EvalInputs:
    """What evaluating a run needs besides its clusters and match; fixed for a claim."""

    spec: ClaimSpec
    counts: Any
    lognorm: pd.DataFrame
    sample: pd.Series | None
    condition: pd.Series | None
    thr: Thresholds
    n_cpus: int
    pair: pd.Series | None


def _evaluate_run(
    inp: _EvalInputs, clusters: pd.Series, m: Match
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """The claim's evaluation of one matched run, and the informational doublet check."""
    ev = _evaluate(
        inp.spec,
        m,
        clusters,
        inp.counts,
        inp.lognorm,
        inp.sample,
        inp.condition,
        inp.thr,
        inp.n_cpus,
        inp.pair,
    )
    doublet_info = None
    if inp.spec.claim_type == "population_exists" and m.found and m.clusters:
        # informational: no cutoff separates doublet mixtures from intermediate populations
        # with margin (bench/RESULTS.md), so this never changes a verdict
        dc = doublet_cluster_check(
            inp.counts,
            clusters,
            m.clusters,
            residual_ratio_max=inp.thr.doublet_residual_ratio_max,
            libsize_ratio_min=inp.thr.doublet_libsize_ratio_min,
        )
        doublet_info = dc.to_dict()
    return ev, doublet_info


def _run_evaluation(
    i: int,
    params: PipelineParams,
    base_params: PipelineParams,
    rr: RunResult,
    m: Match,
    ev: dict[str, Any],
    doublet_info: dict[str, Any] | None,
    runtime_s: float,
) -> RunEvaluation:
    match_info = m.to_dict()
    match_info.update(_doublet_summary(rr.doublet_scores, m.cells))
    match_info["doublet_cluster"] = doublet_info
    effective: dict[str, Any] = dict(rr.effective)
    match_info["doublet_population"] = _doublet_population_summary(
        effective.get("doublet_clusters"), m.clusters if m.found else []
    )
    label = params.label(base_params)
    capped = effective.get("capped") or {}
    if capped:
        label += " [" + ", ".join(f"{k} capped to {v}" for k, v in capped.items()) + "]"
    return RunEvaluation(
        index=i,
        label=label,
        params=params.model_dump(),
        is_baseline=(i == 0),
        is_one_at_a_time=(len(params.diff(base_params)) == 1),
        n_cells=rr.n_cells,
        n_removed_qc=rr.n_removed_qc,
        n_removed_doublet=rr.n_removed_doublet,
        n_clusters=rr.n_clusters,
        runtime_s=round(runtime_s, 3),
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
        effective_params=effective,
    )


# ------------------------------------------------------------------ worker processes
# BLAS and numba are limited to one thread in each worker, so N workers use about N cores.
# OpenMP is left alone: scikit-learn's exact kNN (used by scrublet and the neighbour graph)
# picks its chunking strategy from the OpenMP thread count, and with one thread it can return
# different neighbours for near-ties, which changes doublet calls. Workers keep the parent's
# OpenMP setting so their results match a sequential run.
_THREAD_VARS = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMBA_NUM_THREADS",
)

_EvalOutput = tuple[dict[str, Any], dict[str, Any] | None, float]


@dataclass
class _WorkerState:
    pipe: Pipeline
    inputs: _EvalInputs


_worker: _WorkerState | None = None


def _init_worker(
    counts: Any,
    grid: GridConfig,
    sample: pd.Series | None,
    batch: pd.Series | None,
    mito_prefix: str | None,
    spec: ClaimSpec,
    lognorm: pd.DataFrame,
    condition: pd.Series | None,
    thr: Thresholds,
    pair: pd.Series | None,
) -> None:
    global _worker
    pipe = Pipeline(counts, grid, sample=sample, batch=batch, mito_prefix=mito_prefix)
    # one process per worker already; pydeseq2's own result does not depend on n_cpus
    inputs = _EvalInputs(spec, counts, lognorm, sample, condition, thr, 1, pair)
    _worker = _WorkerState(pipe, inputs)


def _pipeline_task(params: PipelineParams) -> RunResult:
    assert _worker is not None
    return _worker.pipe.run(params)


def _evaluate_task(clusters: pd.Series, m: Match) -> _EvalOutput:
    assert _worker is not None
    t0 = time.perf_counter()
    ev, doublet_info = _evaluate_run(_worker.inputs, clusters, m)
    return ev, doublet_info, time.perf_counter() - t0


@contextmanager
def _single_threaded_blas() -> Iterator[None]:
    """BLAS and numba limited to one thread in this process, as in the worker processes.

    OpenBLAS (numpy's BLAS on Linux) can give slightly different floating-point results for
    different thread counts, so the sequential grid uses one thread like the workers do.
    """
    import numba
    from threadpoolctl import threadpool_limits

    previous = numba.get_num_threads()
    numba.set_num_threads(1)
    try:
        with threadpool_limits(limits=1, user_api="blas"):
            yield
    finally:
        numba.set_num_threads(previous)


@contextmanager
def _worker_pool(workers: int, initargs: tuple[Any, ...]) -> Iterator[ProcessPoolExecutor]:
    """A pool of freshly spawned workers with single-threaded BLAS and numba.

    The thread limits go into the environment the workers start with (the libraries read them
    at import) and the parent's environment is restored when the pool closes. Spawn rather
    than fork: forking a process whose thread pools are already running is unsafe.
    """
    saved = {k: os.environ.get(k) for k in _THREAD_VARS}
    os.environ.update({k: "1" for k in _THREAD_VARS})
    try:
        pool = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_init_worker,
            initargs=initargs,
        )
        try:
            yield pool
        except BaseException:
            pool.shutdown(wait=True, cancel_futures=True)
            raise
        pool.shutdown(wait=True)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _evaluate(
    spec: ClaimSpec,
    m: Match,
    clusters: pd.Series,
    counts: Any,
    lognorm: pd.DataFrame,
    sample: pd.Series | None,
    condition: pd.Series | None,
    thr: Thresholds,
    n_cpus: int,
    pair: pd.Series | None,
) -> dict[str, Any]:
    if spec.claim_type == "population_exists":
        return evaluate_population(m, thr)
    assert spec.comparison is not None and sample is not None and condition is not None
    if spec.claim_type == "composition_shift":
        return evaluate_composition(
            m.cells, clusters.index, sample, condition, spec.comparison, thr, pair=pair
        )
    if not m.found:
        return {
            "holds": False,
            "holds_cell_level": False,
            "replicate_ok": False,
            "evaluable": False,
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
        n_cpus=n_cpus,
        pair=pair,
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
    pair: pd.Series | None,
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
        return evaluate_composition(
            cells, all_cells, sample, condition, spec.comparison, thr, pair=pair
        )
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
        pair=pair,
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


def _doublet_population_summary(
    info: dict[str, Any] | None, clusters: list[str]
) -> dict[str, Any] | None:
    """The doublet-population statistic of the matched clusters, when the step ran in this run.

    Clusters below the test's minimum size have no statistic and are left out; a tested cluster
    with no pair of parent groups maps to None.
    """
    if not info:
        return None
    stats = info.get("statistics") or {}
    return {
        "statistics": {c: stats[c] for c in clusters if c in stats},
        "auroc_max": info.get("auroc_max"),
    }


def _removed_doublet_clusters(e: RunEvaluation) -> tuple[int | None, str]:
    """Cells and clusters removed as doublet populations, for runs.csv.

    Returns (cells, "cluster:statistic;...") or (None, "") when the step did not run.
    """
    info = e.effective_params.get("doublet_clusters")
    if not info:
        return None, ""
    stats = info.get("statistics") or {}
    parts = []
    for c in info.get("removed") or []:
        st = stats.get(c)
        parts.append(f"{c}:{st:.4f}" if st is not None else str(c))
    return int(info.get("n_cells") or 0), ";".join(parts)


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
                "n_removed_qc",
                "n_removed_doublet",
                "n_removed_doublet_cluster",
                "doublet_clusters_removed",
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
            dc_cells, dc_removed = _removed_doublet_clusters(e)
            w.writerow(
                [
                    e.index,
                    e.label,
                    *[e.params[k] for k in param_keys],
                    e.n_cells,
                    e.n_clusters,
                    e.n_removed_qc,
                    e.n_removed_doublet,
                    dc_cells,
                    dc_removed,
                    e.match.get("n_cells"),
                    e.match.get("jaccard"),
                    e.match.get("marker_auroc"),
                    *[e.key_stats.get(k) for k in stat_keys],
                    e.holds,
                    e.holds_cell_level,
                    e.runtime_s,
                ]
            )
