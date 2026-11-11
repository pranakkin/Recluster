"""Turn per-run outcomes into a verdict, with every threshold explicit."""

from __future__ import annotations

from typing import Any

from recluster.config import PipelineParams, Thresholds
from recluster.results import RunEvaluation, Verdict

ROBUST, FRAGILE, NOT_SUPPORTED = "ROBUST", "FRAGILE", "NOT SUPPORTED"


def compute_verdict(
    runs: list[RunEvaluation],
    baseline: PipelineParams,
    thr: Thresholds,
    claim_type: str,
) -> Verdict:
    n = len(runs)
    n_holds = sum(1 for r in runs if r.holds)
    survival = n_holds / n if n else 0.0

    cell_vals = [r.holds_cell_level for r in runs if r.holds_cell_level is not None]
    cell_survival = sum(cell_vals) / len(cell_vals) if cell_vals else None
    rec_vals = [r.recovered for r in runs if r.recovered is not None]
    recovered_survival = sum(rec_vals) / len(rec_vals) if rec_vals else None

    if survival >= thr.robust_min_survival:
        label = ROBUST
    elif survival >= thr.fragile_min_survival:
        label = FRAGILE
    else:
        label = NOT_SUPPORTED

    pseudoreplication = bool(
        claim_type in ("differential_expression", "composition_shift")
        and cell_survival is not None
        and cell_survival >= thr.pseudoreplication_cell_survival
        and survival < thr.pseudoreplication_replicate_survival
    )

    by_choice = _by_choice(runs, baseline)
    breakers = sorted(
        (row for row in by_choice if row["n_runs"] > row["n_holds"]),
        key=lambda r: (-(1 - r["n_holds"] / r["n_runs"]), -(r["n_runs"] - r["n_holds"]), r["axis"]),
    )
    baseline_run = next((r for r in runs if r.is_baseline), None)
    baseline_holds = bool(baseline_run.holds) if baseline_run else False
    oat_flips = [
        {"label": r.label, "holds": r.holds, "params_changed": _changed(r.params, baseline)}
        for r in runs
        if r.is_one_at_a_time and r.holds != baseline_holds
    ]

    summary = _summary(label, n_holds, n, baseline_holds, oat_flips, pseudoreplication, claim_type)
    return Verdict(
        label=label,
        survival=survival,
        n_runs=n,
        n_holds=n_holds,
        baseline_holds=baseline_holds,
        cell_level_survival=cell_survival,
        recovered_survival=recovered_survival,
        pseudoreplication_flag=pseudoreplication,
        by_choice=by_choice,
        breakers=breakers,
        oat_flips=oat_flips,
        summary=summary,
    )


def _changed(params: dict[str, Any], baseline: PipelineParams) -> dict[str, Any]:
    return {k: v for k, v in params.items() if v != getattr(baseline, k)}


def _by_choice(runs: list[RunEvaluation], baseline: PipelineParams) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for axis in PipelineParams.model_fields:
        values = sorted({r.params[axis] for r in runs}, key=lambda v: (str(type(v)), v))
        if len(values) < 2:
            continue
        for v in values:
            subset = [r for r in runs if r.params[axis] == v]
            rows.append(
                {
                    "axis": axis,
                    "value": v,
                    "n_runs": len(subset),
                    "n_holds": sum(1 for r in subset if r.holds),
                    "is_baseline": v == getattr(baseline, axis),
                }
            )
    return rows


def _summary(
    label: str,
    n_holds: int,
    n: int,
    baseline_holds: bool,
    oat_flips: list[dict[str, Any]],
    pseudoreplication: bool,
    claim_type: str,
) -> str:
    parts = [
        f"The claim held in {n_holds} of {n} runs ({100 * n_holds / n:.0f}%)." if n else "No runs."
    ]
    parts.append(
        "It held under the baseline pipeline."
        if baseline_holds
        else "It did not hold under the baseline pipeline."
    )
    if oat_flips:
        flips = "; ".join(
            ", ".join(f"{k}={v}" for k, v in f["params_changed"].items()) for f in oat_flips[:6]
        )
        verb = "broke" if baseline_holds else "rescued"
        parts.append(f"Single changes that {verb} it: {flips}.")
    if pseudoreplication:
        what = "genes" if claim_type == "differential_expression" else "proportions"
        parts.append(
            f"Cell-level statistics call the {what} significant in most runs while replicate-aware "
            "statistics do not: this pattern is characteristic of pseudoreplication."
        )
    return " ".join(parts)
