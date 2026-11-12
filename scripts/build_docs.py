# ruff: noqa: E501  (CSS and HTML template strings)
"""Copy the benchmark reports into docs/ and write a static docs/index.html.

GitHub Pages serves the committed ``docs/`` folder as-is (a ``.nojekyll`` file disables Jekyll).
Nothing here runs on Pages; this script runs locally after ``recluster bench`` and its output is
committed. Usage: ``python scripts/build_docs.py``.
"""

from __future__ import annotations

import json
import shutil
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
REPO = "https://github.com/pranakkin/Recluster"

CSS = """
:root { --ink:#1b1b1b; --muted:#5f5e58; --rule:#d8d6cf; --bg:#fbfaf7; --panel:#fff; --soft:#f3f1ea;
        --accent:#1f4e79; --pass:#2f6b3a; --fail:#9b3b2e; --warn:#8a6a1e; }
body { margin:0; background:var(--bg); color:var(--ink);
       font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
main { max-width: 960px; margin: 0 auto; padding: 36px 24px 72px; }
.eyebrow { font-size: 12px; letter-spacing: .12em; text-transform: uppercase; color: var(--muted); }
h1 { font-size: 27px; font-weight: 650; margin: 6px 0 4px; }
h2 { font-size: 17px; font-weight: 650; margin: 36px 0 10px; padding-bottom: 6px; border-bottom: 1px solid var(--rule); }
p { margin: 8px 0; } a { color: var(--accent); }
table { border-collapse: collapse; width: 100%; font-size: 14px; background: var(--panel); }
th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--rule); vertical-align: top; }
th { font-weight: 600; color: var(--muted); font-size: 12px; text-transform: uppercase; letter-spacing: .05em; background: var(--soft); }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: .92em; white-space: nowrap; }
.v { font-weight: 700; } .v.robust { color: var(--pass); } .v.fragile { color: var(--warn); } .v.not-supported { color: var(--fail); }
.small { font-size: 13px; color: var(--muted); }
"""

ICON = {"ROBUST": "✓", "FRAGILE": "◐", "NOT SUPPORTED": "✕"}


def outcome_text(c: dict) -> str:
    if c["role"] == "control":
        return "control OK" if c["passes"] else "control failed"
    if c["part"] == "planted" or c["role"] == "planted_time":
        return "caught" if c["passes"] else "missed"
    if c["role"] in ("false", "null"):
        return "correct" if c["passes"] else "false alarm"
    return "as expected" if c["passes"] else "missed"


def main() -> None:
    results = json.loads((ROOT / "bench" / "results.json").read_text())
    reports = DOCS / "reports"
    shutil.rmtree(reports, ignore_errors=True)
    rows: list[tuple[str, dict, str]] = []
    for c in results["cases"]:
        src = Path(c["report"])
        if not src.is_absolute():
            src = ROOT / src
        if not src.exists():
            raise SystemExit(f"missing report {src}; run `recluster bench` first")
        rel = f"reports/{c['part']}/{src.parent.name}.html"
        dst = DOCS / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        rows.append((c["part"], c, rel))

    def table(part: str, title: str, blurb: str) -> str:
        body = []
        for i, (p, c, rel) in enumerate(rows, 1):
            if p != part:
                continue
            v = c["verdict"]
            body.append(
                f"<tr><td class='mono'>{i}</td><td><a href='{rel}'>{escape(c['title'])}</a>"
                f"<div class='small'>{escape(c['claim'])}</div></td>"
                f"<td class='small'>{escape(c['expected'])}</td>"
                f"<td class='v {v.lower().replace(' ', '-')}'>{ICON[v]} {v}</td>"
                f"<td class='mono'>{c['n_holds']}/{c['n_runs']}</td><td>{outcome_text(c)}</td></tr>"
            )
        return (
            f"<h2>{escape(title)}</h2><p class='small'>{blurb}</p>"
            "<table><thead><tr><th>#</th><th>Case</th><th>Expected</th><th>Verdict</th>"
            f"<th>Held</th><th>Outcome</th></tr></thead><tbody>{''.join(body)}</tbody></table>"
        )

    n_cases = len(rows)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Recluster · benchmark reports</title><style>{CSS}</style></head>
<body><main>
<div class="eyebrow">Recluster · benchmark reports</div>
<h1>An independent stress test for single-cell RNA-seq claims</h1>
<p>{
        n_cases
    } claims from the benchmark, each with its full report: verdict, survival across the grid,
what broke it, replicate-aware vs cell-level statistics, methods and provenance. Numbers are from
the run of {escape(results["generated_at"][:10])}; <a href="{
        REPO
    }/blob/main/bench/RESULTS.md">bench/RESULTS.md</a>
explains each miss and the false alarm; <a href="{
        REPO
    }/blob/main/docs/HOW_IT_WORKS.md">HOW_IT_WORKS.md</a>
describes the design. Source: <a href="{REPO}">{REPO}</a>.</p>
{
        table(
            "kang",
            "Real data: Kang et al. 2018",
            "PBMCs from 8 lupus patients, IFN-β vs control. True claims, a claim with no known effect, "
            "claims that are false by construction (fake donor groups, donor-level label shuffles), "
            "and a planted time-of-day effect claimed as a condition effect.",
        )
    }
{
        table(
            "planted",
            "Planted errors on PBMC3k",
            "Each error is built deterministically from scanpy's PBMC3k and claimed as a naive analysis would; "
            "the caught/missed rule was fixed before running.",
        )
    }
<p class="small">Reports are self-contained HTML produced by <span class="mono">recluster bench</span>; this page is a static
index with no build step.</p>
</main></body></html>
"""
    (DOCS / "index.html").write_text(html)
    (DOCS / ".nojekyll").write_text("")
    print(f"wrote docs/index.html and {n_cases} reports under docs/reports/")


if __name__ == "__main__":
    main()
