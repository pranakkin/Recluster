"""Score the ground-truth claims against the bulk truth, as PLAN.md "Scoring" specifies.

Four callers, each scored on the same claims:
- Recluster ROBUST (the call), and ROBUST or FRAGILE (also reported);
- naive cell-level testing: the baseline run's cell-level test (Wilcoxon, `holds_cell_level`);
- a single pseudobulk analysis: the baseline run's replicate-level test alone (`holds`).

Sensitivity = calls among true-change claims, specificity = no call among no-change claims, each
with a Wilson 95 % interval, per comparison and pooled. Writes out/scores.json and
out/scores.md and lists every miss and false call.

Usage: python bench/ground_truth/score.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "out"
CALLERS = ("robust", "robust_or_fragile", "naive_cell", "single_pseudobulk")
NAMES = {
    "robust": "Recluster ROBUST",
    "robust_or_fragile": "Recluster ROBUST or FRAGILE",
    "naive_cell": "naive cell-level test (baseline run)",
    "single_pseudobulk": "one pseudobulk analysis (baseline run)",
}


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, centre - half), min(1.0, centre + half))


def calls_for(prov: dict) -> dict[str, bool]:
    label = prov["verdict"]["label"]
    base = next(r for r in prov["runs"] if r["is_baseline"])
    return {
        "robust": label == "ROBUST",
        "robust_or_fragile": label in ("ROBUST", "FRAGILE"),
        "naive_cell": bool(base.get("holds_cell_level")),
        "single_pseudobulk": bool(base.get("holds")),
    }


def rate(rows: list[dict], cls: str, caller: str) -> dict:
    sub = [r for r in rows if r["truth"] == cls]
    if cls == "change":  # sensitivity: called
        k = sum(r["calls"][caller] for r in sub)
    else:  # specificity: not called
        k = sum(not r["calls"][caller] for r in sub)
    lo, hi = wilson(k, len(sub))
    return {"k": k, "n": len(sub), "rate": k / len(sub) if sub else None, "ci95": [lo, hi]}


def main() -> None:
    claims = json.loads((HERE / "claims" / "claims.json").read_text())
    summary = json.loads((OUT / "summary.json").read_text())
    rows = []
    for c in claims:
        key = f"{c['comparison']}/{c['gene']}"
        s = summary.get(key)
        if s is None or "error" in s:
            continue
        prov = json.loads((HERE / s["folder"] / "provenance.json").read_text())
        base = next(r for r in prov["runs"] if r["is_baseline"])
        rows.append(
            {
                "key": key,
                "comparison": c["comparison"],
                "gene": c["gene"],
                "truth": c["truth"],
                "bulk_log2fc": c["bulk_log2fc"],
                "bulk_padj": c["bulk_padj"],
                "verdict": prov["verdict"]["label"],
                "n_holds": prov["verdict"]["n_holds"],
                "n_runs": prov["verdict"]["n_runs"],
                "baseline_stats": base.get("key_stats"),
                "calls": calls_for(prov),
            }
        )
    groups = {"pooled": rows}
    for comp in sorted({r["comparison"] for r in rows}):
        groups[comp] = [r for r in rows if r["comparison"] == comp]
    scores = {
        g: {
            caller: {
                "sensitivity": rate(rs, "change", caller),
                "specificity": rate(rs, "no_change", caller),
            }
            for caller in CALLERS
        }
        for g, rs in groups.items()
    }
    errors = {
        caller: {
            "misses": [r["key"] for r in rows if r["truth"] == "change" and not r["calls"][caller]],
            "false_calls": [
                r["key"] for r in rows if r["truth"] == "no_change" and r["calls"][caller]
            ],
        }
        for caller in CALLERS
    }
    result = {"n_scored": len(rows), "scores": scores, "errors": errors, "claims": rows}
    (OUT / "scores.json").write_text(json.dumps(result, indent=1))
    md = tables(result)
    (OUT / "scores.md").write_text(md)
    print(md)


def fmt(x: dict) -> str:
    if not x["n"]:
        return "0/0"
    lo, hi = x["ci95"]
    return f"{x['k']}/{x['n']} ({100 * x['rate']:.0f} %, {100 * lo:.0f}-{100 * hi:.0f})"


def tables(result: dict) -> str:
    lines = [f"Scored claims: {result['n_scored']}", ""]
    for g, by_caller in result["scores"].items():
        lines += [f"### {g}", "", "| caller | sensitivity (95 % CI) | specificity (95 % CI) |"]
        lines.append("|---|---|---|")
        for caller, s in by_caller.items():
            lines.append(f"| {NAMES[caller]} | {fmt(s['sensitivity'])} | {fmt(s['specificity'])} |")
        lines.append("")
    lines += ["### Misses and false calls", ""]
    for caller, e in result["errors"].items():
        lines.append(f"- {NAMES[caller]}: misses {e['misses'] or 'none'}; false calls ")
        lines[-1] += f"{e['false_calls'] or 'none'}"
    lines += ["", "### Every claim", "", "| claim | truth | bulk log2FC | verdict | runs | "
              "naive cell | one pseudobulk |", "|---|---|---|---|---|---|---|"]  # fmt: skip
    for r in result["claims"]:
        c = r["calls"]
        lines.append(
            f"| {r['key']} | {r['truth']} | {r['bulk_log2fc']:.2f} | {r['verdict']} | "
            f"{r['n_holds']}/{r['n_runs']} | {'call' if c['naive_cell'] else '-'} | "
            f"{'call' if c['single_pseudobulk'] else '-'} |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
