"""Collect the audit outputs into out/results.json and out/results_tables.md.

Reports, for every spec, Recluster's verdict; for DE specs, each gene's own survival from the
same runs (PROTOCOL amendment 8); for the 11 abundance specs, the per-run BH adjustment over the
11 types (amendment 7). The last two are our analysis of Recluster's per-run output, not
Recluster verdicts.

Usage: uv run python bench/agent_audit/summarize.py
"""

from __future__ import annotations

import json
import statistics
from pathlib import Path

import numpy as np

from recluster.evaluators.de import bh_adjust

HERE = Path(__file__).parent
OUT = HERE / "out"
ROBUST_MIN = 0.80


def per_gene_survival(prov: dict) -> dict[str, int]:
    """Runs in which each claimed gene holds at the replicate level (unmatched runs count as no)."""
    genes = prov["claim"]["genes"]
    counts = dict.fromkeys(genes, 0)
    for run in prov["runs"]:
        details = run.get("details") or {}
        if not details.get("evaluable"):
            continue
        for g in details.get("per_gene", []):
            if g.get("holds_replicate"):
                counts[g["gene"]] += 1
    return counts


def composition_bh(provs: dict[str, dict]) -> dict:
    """BH over the abundance specs within each run; a missing p (unmatched) counts as 1."""
    labels = sorted(provs)
    n_runs = {len(p["runs"]) for p in provs.values()}
    assert len(n_runs) == 1, "abundance specs have different grids"
    (n,) = n_runs
    passes = dict.fromkeys(labels, 0)
    min_padj: dict[str, list[float]] = {lab: [] for lab in labels}
    for i in range(n):
        ps = []
        for lab in labels:
            run = provs[lab]["runs"][i]
            assert run["index"] == i
            p = (run.get("key_stats") or {}).get("p_replicate")
            ps.append(1.0 if p is None or not (run.get("details") or {}).get("evaluable") else p)
        padj = bh_adjust(np.asarray(ps, dtype=float))
        for lab, q in zip(labels, padj, strict=True):
            min_padj[lab].append(float(q))
            passes[lab] += int(q < 0.05)
    return {
        "n_runs": n,
        "passes_after_bh": passes,
        "median_padj": {lab: statistics.median(v) for lab, v in min_padj.items()},
        "contradicted": [lab for lab in labels if passes[lab] >= ROBUST_MIN * n],
    }


def main() -> None:
    items = json.loads((HERE / "claims" / "claims.json").read_text())
    summary = json.loads((OUT / "summary.json").read_text())
    rows, comp_provs = [], {}
    for c in items:
        s = summary.get(c["id"])
        row = {k: c[k] for k in ("id", "statement", "agent_says", "secondary")}
        row["superseded_by"] = next(
            (d["id"] for d in items if d.get("supersedes") == c["id"]), None
        )  # amendment 9
        row["claim_type"] = c["spec"]["claim_type"]
        if s is None or "error" in s:
            row["verdict"] = "NOT RUN" if s is None else "ERROR"
            row["error"] = (s or {}).get("error")
            rows.append(row)
            continue
        prov = json.loads((HERE / s["folder"] / "provenance.json").read_text())
        v = prov["verdict"]
        row.update(
            verdict=v["label"],
            label_before_time_check=v.get("label_before_time_check"),
            time_check=v.get("time_check_label"),
            n_holds=v["n_holds"],
            n_runs=v["n_runs"],
            n_not_testable=v.get("n_not_testable"),
            cell_level_survival=v.get("cell_level_survival"),
            breakers=v.get("breakers"),
            report=f"{s['folder']}/report.html",
        )
        jac = [r["match"].get("jaccard") for r in prov["runs"] if r.get("match")]
        jac = [j for j in jac if j is not None]
        row["median_jaccard"] = statistics.median(jac) if jac else None
        if row["claim_type"] == "differential_expression":
            row["per_gene_holds"] = per_gene_survival(prov)
        if row["claim_type"] == "composition_shift":
            comp_provs[c["spec"]["target"]["label_value"]] = prov
            base = next(r for r in prov["runs"] if r["is_baseline"])
            d = base.get("details") or {}
            row["baseline"] = {
                k: d.get(k)
                for k in (
                    "mean_proportion_a",
                    "mean_proportion_b",
                    "p_logit_ttest",
                    "p_mannwhitney",
                )
            }
        rows.append(row)
    result = {"claims": rows}
    if len(comp_provs) == 11:
        result["abundance_bh"] = composition_bh(comp_provs)
    (OUT / "results.json").write_text(json.dumps(result, indent=1, default=str))
    (OUT / "results_tables.md").write_text(tables(result))
    print(tables(result))


def tables(result: dict) -> str:
    lines = ["| id | agent says | verdict | runs | per-gene survival (genes below all runs) |"]
    lines.append("|---|---|---|---|---|")
    for r in result["claims"]:
        if r["verdict"] in ("NOT RUN", "ERROR"):
            lines.append(f"| {r['id']} | {r['agent_says']} | {r['verdict']} | | |")
            continue
        genes = r.get("per_gene_holds") or {}
        weak = ", ".join(f"{g} {k}/{r['n_runs']}" for g, k in genes.items() if k < r["n_runs"])
        tag = " (secondary)" if r["secondary"] else ""
        if r.get("superseded_by"):
            tag = f" (first run, superseded by {r['superseded_by']})"
        gene_col = (weak or f"all genes {r['n_runs']}") if genes else ""
        lines.append(
            f"| {r['id']}{tag} | {r['agent_says']} | {r['verdict']} | "
            f"{r['n_holds']}/{r['n_runs']} | {gene_col} |"
        )
    bh = result.get("abundance_bh")
    if bh:
        lines += ["", "Abundance, BH over the 11 types within each run (runs passing / total):", ""]
        lines += [f"- {lab}: {k}/{bh['n_runs']}" for lab, k in bh["passes_after_bh"].items()]
        lines.append(f"- contradicted (>= 80 % of runs): {bh['contradicted'] or 'none'}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
