# ruff: noqa: E501  (CSS and HTML template strings)
"""Copy the benchmark reports into docs/ and write a static docs/benchmark.html.

GitHub Pages serves the committed ``docs/`` folder as-is (a ``.nojekyll`` file disables Jekyll).
Nothing here runs on Pages; this script runs locally after ``recluster bench`` and its output is
committed. Usage: ``python scripts/build_docs.py``.
"""

from __future__ import annotations

import json
import re
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
ul { padding-left: 20px; } li { margin: 4px 0; }
@media (max-width: 700px) { .wide table { display: block; overflow-x: auto; } }
"""

ICON = {"ROBUST": "✓", "FRAGILE": "◐", "NOT SUPPORTED": "✕"}

# The agent audit reports are not part of results.json. Their published copies under
# docs/reports/agent/ had local paths removed by hand, so they are never re-copied from
# bench/agent_audit/out; this script only adds the navigation banner to them in place.
AGENT_REPORTS = (
    (
        "agent-megakaryocyte-platelet-is-a-distinct-population-per-sa.html",
        "Megakaryocyte/platelet population",
        "NOT SUPPORTED",
        "1/21",
    ),
    ("agent-nk-down-in-nk-nell2.html", "NELL2 down in NK", "FRAGILE", "6/21"),
    ("agent-myeloid-down-in-cdc.html", "down-genes in cDC", "FRAGILE", "5/21"),
    (
        "agent-null-cd14-mono-proportion-differs-stim-vs-ctrl.html",
        "CD14 monocyte abundance",
        "FRAGILE",
        "15/21",
    ),
)

NAV_START, NAV_END = "<!-- site-nav -->", "<!-- /site-nav -->"
CSS_START, CSS_END = "/* site-nav */", "/* /site-nav */"
NAV_CSS = (
    ".site-nav a{color:#1f4e79}@media (max-width:700px){table{display:block;overflow-x:auto}"
    "dd.mono{white-space:normal;overflow-wrap:anywhere}}"
)
NAV_LINE = (
    '<p class="site-nav" style="font:13px/1.5 sans-serif;margin:0 0 16px"><a href="../../index.html">Recluster</a>'
    ' · <a href="../../benchmark.html">All benchmark reports</a>'
    f' · <a href="{REPO}/blob/main/docs/HOW_IT_WORKS.md#how-to-read-a-report">How to read a report</a></p>'
)


def outcome_text(c: dict) -> str:
    if c["role"] == "control":
        return "control OK" if c["passes"] else "control failed"
    if c["part"] == "planted" or c["role"] == "planted_time":
        return "caught" if c["passes"] else "missed"
    if c["role"] in ("false", "null"):
        return "correct" if c["passes"] else "false alarm"
    return "as expected" if c["passes"] else "missed"


def add_banner(path: Path, context: str) -> None:
    """Insert (or replace) the back-link banner after <main> and the narrow-screen table CSS.

    Idempotent: an earlier banner between the site-nav markers is removed first.
    """
    html = path.read_text()
    html = re.sub(re.escape(NAV_START) + ".*?" + re.escape(NAV_END) + "\n?", "", html, flags=re.S)
    html = re.sub(re.escape(CSS_START) + ".*?" + re.escape(CSS_END) + "\n?", "", html, flags=re.S)
    banner = f"{NAV_START}{NAV_LINE}"
    if context:
        banner += f'<p class="site-nav" style="font:13px/1.5 sans-serif;margin:-10px 0 16px;color:#5f5e58">{context}</p>'
    banner += NAV_END + "\n"
    if "<main>" not in html or "</style>" not in html:
        raise SystemExit(f"{path}: no <main> or </style> to anchor the banner")
    html = html.replace("<main>\n", "<main>\n" + banner, 1)
    html = html.replace("</style>", f"{CSS_START}{NAV_CSS}{CSS_END}\n</style>", 1)
    path.write_text(html)


def main() -> None:
    results = json.loads((ROOT / "bench" / "results.json").read_text())
    reports = DOCS / "reports"
    # Remove only the parts this script rebuilds; reports/agent/ is maintained by hand.
    for part in sorted({c["part"] for c in results["cases"]}):
        shutil.rmtree(reports / part, ignore_errors=True)
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

    n_total = len(rows)
    for i, (_, c, rel) in enumerate(rows, 1):
        context = (
            f"Benchmark case {i} of {n_total} · Expected: {escape(c['expected'])} · "
            f"Recluster: {escape(c['verdict'])}, {c['n_holds']}/{c['n_runs']} runs · "
            f"{outcome_text(c)}"
        )
        add_banner(DOCS / rel, context)
    agent_details = f"{REPO}/blob/main/bench/agent_audit/RESULTS.md"
    for name, *_ in AGENT_REPORTS:
        path = reports / "agent" / name
        if not path.exists():
            raise SystemExit(f"missing {path}; restore docs/reports/agent/ from git")
        add_banner(
            path, f'Agent audit case study (Kang 2018). <a href="{agent_details}">Details</a>'
        )

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
    agent_items = "".join(
        f"<li><a href='reports/agent/{name}'>{escape(text)}</a> "
        f"<span class='v {v.lower().replace(' ', '-')}'>{ICON[v]} {v}</span>, "
        f"<span class='mono'>{held}</span></li>"
        for name, text, v, held in AGENT_REPORTS
    )
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Recluster benchmark: {n_cases} claims with known answers</title><style>{CSS}</style></head>
<body><main>
<div class="eyebrow"><a href="./">← Project page</a></div>
<h1>Recluster benchmark: {n_cases} claims with known answers</h1>
<p>Each row is one claim with a known answer. <b>Expected</b>: the result the pass rule asks for,
fixed before the run. <b>Verdict</b>: what Recluster returned (✓ ROBUST, ◐ FRAGILE, ✕ NOT SUPPORTED).
<b>Held</b>: re-analyses in which the claim held, out of the grid. <b>Outcome</b>: whether that
passes the rule (as expected / correct / caught / control OK = pass; false alarm / missed = fail).
Click a case for its full report.</p>
<p>Numbers are from the run of {escape(results["generated_at"][:10])}; <a href="{
        REPO
    }/blob/main/bench/RESULTS.md">bench/RESULTS.md</a>
explains each miss and the false alarm; <a href="{
        REPO
    }/blob/main/docs/HOW_IT_WORKS.md">HOW_IT_WORKS.md</a>
describes the design and <a href="{
        REPO
    }/blob/main/docs/HOW_IT_WORKS.md#how-to-read-a-report">how to read a report</a>. Source: <a href="{
        REPO
    }">{REPO}</a>.</p>
<div class="wide">
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
</div>
<h2 id="agent">Agent audit: Kang 2018</h2>
<p class='small'>An analysis agent's full analysis of Kang 2018, checked claim by claim. Of its 36 primary
population and expression claims, 32 came out ROBUST, 3 FRAGILE and 1 NOT SUPPORTED. The
population results come from a re-run with per-donor doublet removal, decided after the first
results had been seen (both runs are reported). The agent found no change in abundance; the
FRAGILE abundance report below does not contradict that once the 11 types are corrected for
together.
<a href="{REPO}/blob/main/bench/agent_audit/RESULTS.md">bench/agent_audit/RESULTS.md</a> has the
protocol and every result. Four of the reports:</p>
<ul>{agent_items}</ul>
<p class="small">Reports are self-contained HTML produced by <span class="mono">recluster bench</span> and
<span class="mono">recluster audit</span>; this page is a static index with no build step.</p>
</main></body></html>
"""
    (DOCS / "benchmark.html").write_text(html)
    (DOCS / ".nojekyll").write_text("")
    print(
        f"wrote docs/benchmark.html, {n_cases} reports under docs/reports/ and the banner on "
        f"{len(AGENT_REPORTS)} agent reports"
    )


if __name__ == "__main__":
    main()
