"""Run the whole benchmark and write machine-readable results.

``recluster bench`` (or ``python -m recluster.bench.runner``) builds every case, runs each claim
through the normal engine with the default grid and thresholds, applies each case's
pre-registered pass/caught rule, and writes ``bench/results.json`` plus a Markdown table.
Nothing in here looks at a result before deciding what to run.
"""

from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from recluster.bench.doublets import doublet_check
from recluster.bench.kang import (
    check_kang,
    donor_level_shuffle,
    download_kang,
    fake_split_within_control,
    kang_cases,
    load_kang,
)
from recluster.bench.pbmc3k import load_pbmc3k, reference_labels
from recluster.bench.planted import (
    plant_batch_effect,
    plant_doublets,
    plant_overclustering,
    plant_qc_deletion,
    scrublet_recall_on_injected,
)
from recluster.config import GridConfig, Thresholds
from recluster.engine import run_claim, slugify
from recluster.spec import ClaimSpec

BENCH_DIR = Path(__file__).resolve().parents[2] / "bench"


def _say(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def _run_one(job: dict[str, Any]) -> dict[str, Any]:
    """Executed in a worker process: run one claim, return a compact summary."""
    spec = ClaimSpec.model_validate(job["spec"])
    grid = GridConfig.model_validate(job["grid"])
    thr = Thresholds.model_validate(job["thresholds"])
    t0 = time.perf_counter()
    res = run_claim(
        job["data_path"],
        spec,
        out_dir=job["out_dir"],
        grid=grid,
        thresholds=thr,
        audit=job["audit"],
        progress=None,
    )
    v = res.verdict
    holding = [r for r in res.runs if r.holds]
    return {
        "key": job["key"],
        "verdict": v.label,
        "n_holds": v.n_holds,
        "n_runs": v.n_runs,
        "survival": round(v.survival, 4),
        "cell_level_survival": None
        if v.cell_level_survival is None
        else round(v.cell_level_survival, 4),
        "pseudoreplication_flag": v.pseudoreplication_flag,
        "baseline_holds": v.baseline_holds,
        "top_breakers": [
            {
                "axis": b["axis"],
                "value": b["value"],
                "failed": b["n_runs"] - b["n_holds"],
                "runs": b["n_runs"],
            }
            for b in v.breakers[:3]
        ],
        "oat_flips": [f["label"] for f in v.oat_flips],
        "holding_runs_params": [r.params for r in holding],
        "runs": [
            {
                "index": r.index,
                "label": r.label,
                "params": r.params,
                "holds": r.holds,
                "holds_cell_level": r.holds_cell_level,
            }
            for r in res.runs
        ],
        "as_labelled": _compact_as_labelled(res.as_labelled),
        "report": str(Path(job["out_dir"]) / slugify(spec.name) / "report.html"),
        "runtime_s": round(time.perf_counter() - t0, 1),
    }


def _compact_as_labelled(d: dict[str, Any] | None) -> dict[str, Any] | None:
    if not d:
        return None
    keep = (
        "holds",
        "holds_cell_level",
        "p_logit_ttest",
        "p_mannwhitney",
        "p_cell_level_chi2",
        "n_samples_a",
        "n_samples_b",
    )
    out = {k: d[k] for k in keep if k in d}
    if "per_gene" in d:
        out["per_gene"] = [
            {
                k: g.get(k)
                for k in ("gene", "cell_padj", "pb_padj", "cell_log2fc", "pb_log2fc", "agreement")
            }
            for g in d["per_gene"]
        ]
    return out


class _ResultView:
    """Minimal object exposing what the pre-registered rules look at."""

    def __init__(self, summary: dict[str, Any]) -> None:
        from types import SimpleNamespace

        self.verdict = SimpleNamespace(
            label=summary["verdict"],
            breakers=[{"axis": b["axis"]} for b in summary["top_breakers"]],
        )
        self.runs = [SimpleNamespace(holds=r["holds"], params=r["params"]) for r in summary["runs"]]


def build_kang_jobs(
    data_dir: Path, out_dir: Path, grid: GridConfig, thr: Thresholds
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    real = download_kang(data_dir / "kang_2018.h5ad")
    adata = load_kang(real)
    info = check_kang(adata)
    split_path = data_dir / "kang_fake_split_ctrl.h5ad"
    shuffle_path = data_dir / "kang_fake_shuffle.h5ad"
    split = fake_split_within_control(adata, seed=0)
    split.write_h5ad(split_path)
    shuffled = donor_level_shuffle(adata, seed=0)
    shuffled.write_h5ad(shuffle_path)
    cases = kang_cases(real, split_path, shuffle_path, dict(split.uns["fake_split"]))
    jobs = [
        {
            "key": c.key,
            "spec": c.spec.model_dump(exclude_none=True),
            "data_path": str(c.data_path),
            "out_dir": str(out_dir / "kang"),
            "grid": grid.model_dump(),
            "thresholds": thr.model_dump(),
            "audit": False,
        }
        for c in cases
    ]
    meta = {
        c.key: {
            "title": c.title,
            "how_built": c.how_built,
            "expected": c.expected,
            "rule": c.passes_rule,
            "role": c.role,
            "claim": c.spec.name,
            "passes": c.passes,
        }
        for c in cases
    }
    info["fake_split"] = dict(split.uns["fake_split"])
    info["fake_shuffle"] = dict(shuffled.uns["fake_shuffle"])
    _say("doublet check against demuxlet calls (GEO GSE96583 batch 2)")
    info["doublet_check"] = doublet_check(data_dir / "geo")
    return jobs, meta, info


def build_planted_jobs(
    data_dir: Path, out_dir: Path, grid: GridConfig, thr: Thresholds
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    adata = load_pbmc3k(data_dir)
    ref = reference_labels(adata)
    info: dict[str, Any] = {"reference_label_counts": ref.value_counts().to_dict()}
    cases = []
    cases += plant_doublets(adata, ref, seed=0)
    info["scrublet_on_injected"] = scrublet_recall_on_injected(cases[0].adata)
    cases += plant_overclustering(adata, ref, seed=0)
    cases += plant_qc_deletion(adata, ref, seed=0)
    cases += plant_batch_effect(adata, ref, seed=0)
    planted_dir = data_dir / "planted"
    planted_dir.mkdir(parents=True, exist_ok=True)
    jobs, meta = [], {}
    written: dict[int, Path] = {}
    for c in cases:
        if id(c.adata) not in written:
            p = planted_dir / f"pbmc3k_{c.key}.h5ad"
            c.adata.write_h5ad(p)
            written[id(c.adata)] = p
        jobs.append(
            {
                "key": c.key,
                "spec": c.spec.model_dump(exclude_none=True),
                "data_path": str(written[id(c.adata)]),
                "out_dir": str(out_dir / "planted"),
                "grid": grid.model_dump(),
                "thresholds": thr.model_dump(),
                "audit": c.audit,
            }
        )
        meta[c.key] = {
            "title": c.title,
            "how_built": c.how_built,
            "expected": c.expected,
            "rule": c.caught_rule,
            "role": c.role,
            "claim": c.spec.name,
            "passes": c.caught,
            "extra": c.extra,
        }
    return jobs, meta, info


def run_bench(
    *,
    data_dir: Path = Path("data"),
    out_dir: Path = BENCH_DIR / "out",
    parts: tuple[str, ...] = ("kang", "planted"),
    jobs: int = 4,
    grid: GridConfig | None = None,
    thresholds: Thresholds | None = None,
) -> dict[str, Any]:
    grid = grid or GridConfig()
    thr = thresholds or Thresholds()
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    all_jobs: list[dict[str, Any]] = []
    meta: dict[str, Any] = {}
    info: dict[str, Any] = {}
    if "kang" in parts:
        _say("building Kang 2018 cases")
        j, m, i = build_kang_jobs(data_dir, out_dir, grid, thr)
        all_jobs += j
        meta.update(m)
        info["kang"] = i
    if "planted" in parts:
        _say("building planted-error cases on PBMC3k")
        j, m, i = build_planted_jobs(data_dir, out_dir, grid, thr)
        all_jobs += j
        meta.update(m)
        info["planted"] = i
    _say(f"running {len(all_jobs)} claims with {jobs} worker(s)")

    summaries: dict[str, dict[str, Any]] = {}
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(_run_one, job): job["key"] for job in all_jobs}
        for fut in as_completed(futures):
            s = fut.result()
            summaries[s["key"]] = s
            _say(
                f"  {s['key']:<28s} {s['verdict']:<14s} {s['n_holds']:>2}/{s['n_runs']}  "
                f"{s['runtime_s']:.0f}s"
            )

    rows = []
    for job in all_jobs:
        key = job["key"]
        s, m = summaries[key], meta[key]
        passes = bool(m["passes"](_ResultView(s)))
        rows.append(
            {
                "key": key,
                "part": "kang" if key.startswith(("kang", "fake")) else "planted",
                "role": m["role"],
                "title": m["title"],
                "claim": m["claim"],
                "claim_type": job["spec"]["claim_type"],
                "how_built": m["how_built"],
                "expected": m["expected"],
                "rule": m["rule"],
                "passes": passes,
                "extra": m.get("extra", {}),
                **{k: v for k, v in s.items() if k not in ("key",)},
            }
        )
    results = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "total_runtime_s": round(time.perf_counter() - started, 1),
        "grid": grid.model_dump(),
        "thresholds": thr.model_dump(),
        "datasets": info,
        "cases": rows,
    }
    BENCH_DIR.mkdir(exist_ok=True)
    (BENCH_DIR / "results.json").write_text(json.dumps(results, indent=2, default=str))
    (BENCH_DIR / "results_table.md").write_text(results_table(results))
    _say(
        f"wrote {BENCH_DIR / 'results.json'} and results_table.md "
        f"in {results['total_runtime_s']:.0f}s"
    )
    return results


def results_table(results: dict[str, Any]) -> str:
    """Markdown table of outcomes exactly as measured."""
    lines = [
        "| # | Case | Claim | Expected | Verdict | Held | Cell-level | Outcome |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for i, r in enumerate(results["cases"], 1):
        cell = (
            "–"
            if r["cell_level_survival"] is None or r.get("claim_type") == "population_exists"
            else f"{r['cell_level_survival']:.0%}"
        )
        if r["role"] == "control":
            outcome = "control OK" if r["passes"] else "control FAILED"
        elif r["role"] in ("false", "null"):
            outcome = (
                "correct"
                if r["passes"]
                else "FALSE ALARM"
                if r["verdict"] != "NOT SUPPORTED"
                else "correct"
            )
        elif r["part"] == "planted":
            outcome = "caught" if r["passes"] else "MISSED"
        else:
            outcome = "correct" if r["passes"] else "MISSED"
        flag = " ⚑" if r["pseudoreplication_flag"] else ""
        lines.append(
            f"| {i} | {r['title']} | {r['claim']} | {r['expected']} | **{r['verdict']}** | "
            f"{r['n_holds']}/{r['n_runs']} | {cell}{flag} | {outcome} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    import argparse

    p = argparse.ArgumentParser(description="Run the Recluster benchmark.")
    p.add_argument("--data-dir", default="data", type=Path)
    p.add_argument("--out", default=BENCH_DIR / "out", type=Path)
    p.add_argument("--part", choices=["kang", "planted", "all"], default="all")
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args(argv)
    parts = ("kang", "planted") if a.part == "all" else (a.part,)
    run_bench(data_dir=a.data_dir, out_dir=a.out, parts=parts, jobs=a.jobs)


if __name__ == "__main__":  # pragma: no cover
    main()
