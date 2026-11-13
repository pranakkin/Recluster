"""Write docs/site-data.js: the per-run outcomes the project page draws as survival strips.

Reads bench/results.json and each case's provenance.json (written by `recluster bench`), so every
square on the page is a recorded grid run. Usage: ``python scripts/build_site_data.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

from build_docs import outcome_text

ROOT = Path(__file__).resolve().parents[1]
STAT_KEYS = ("pb_log2fc", "pb_padj", "cell_padj", "p_replicate", "jaccard")
THRESHOLDS = ("robust_min_survival", "fragile_min_survival", "alpha")


def run_record(r: dict) -> dict:
    match = r.get("match") or {}
    found = match.get("found", True)
    stats = {k: v for k, v in (r.get("key_stats") or {}).items() if k in STAT_KEYS}
    if match.get("jaccard") is not None:
        stats["jaccard"] = match["jaccard"]
    return {
        "i": r["index"],
        "label": r["label"],
        "state": "held" if r["holds"] else ("nomatch" if not found else "failed"),
        "cell": r.get("holds_cell_level"),
        "cells": r.get("n_cells"),
        "clusters": r.get("n_clusters"),
        "stats": stats,
    }


def main() -> None:
    results = json.loads((ROOT / "bench" / "results.json").read_text())
    cases = []
    for n, c in enumerate(results["cases"], 1):
        report = Path(c["report"])
        prov_path = (report if report.is_absolute() else ROOT / report).parent / "provenance.json"
        prov = json.loads(prov_path.read_text())
        v = prov["verdict"]
        assert (v["label"], v["n_holds"], v["n_runs"]) == (c["verdict"], c["n_holds"], c["n_runs"])
        cases.append(
            {
                "n": n,
                "key": c["key"],
                "part": c["part"],
                "role": c["role"],
                "title": c["title"],
                "claim": c["claim"],
                "expected": c["expected"],
                "verdict": v["label"],
                "before_time_check": v.get("label_before_time_check"),
                "time_check": v.get("time_check_label"),
                "n_holds": v["n_holds"],
                "n_runs": v["n_runs"],
                "outcome": outcome_text(c),
                "report": f"reports/{c['part']}/{report.parent.name}.html",
                "runs": [run_record(r) for r in sorted(prov["runs"], key=lambda r: r["index"])],
            }
        )
    data = {
        "generated_at": results["generated_at"],
        "thresholds": {k: results["thresholds"][k] for k in THRESHOLDS},
        "axes": results["grid"]["axes"],
        "cases": cases,
    }
    out = ROOT / "docs" / "site-data.js"
    out.write_text("window.RECLUSTER = " + json.dumps(data, separators=(",", ":")) + ";\n")
    print(f"wrote {out.relative_to(ROOT)} ({len(cases)} cases, {out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
