"""Render a ClaimResult into one self-contained HTML file."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

from recluster.provenance import json_default, round_floats
from recluster.report.charts import run_strip_svg, survival_by_choice_svg
from recluster.results import ClaimResult

TEMPLATES = Path(__file__).parent / "templates"
ICONS = {"ROBUST": "✓", "FRAGILE": "◐", "NOT SUPPORTED": "✕"}


def _env() -> Environment:
    env = Environment(
        loader=FileSystemLoader(str(TEMPLATES)),
        autoescape=select_autoescape(["html", "j2"]),
        trim_blocks=False,
        lstrip_blocks=False,
    )
    env.filters["thousands"] = lambda v: f"{int(v):,}" if v is not None else "–"
    env.filters["fmtnum"] = fmtnum
    env.filters["fmtval"] = fmtval
    env.filters["params"] = lambda d: ", ".join(f"{k}={fmtval(v)}" for k, v in d.items())
    return env


def fmtnum(v: Any) -> str:
    if v is None:
        return "–"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f != f:
        return "–"
    if abs(f) < 1e-3 and f != 0:
        return f"{f:.1e}"
    return f"{f:.3f}"


def fmtval(v: Any) -> str:
    if isinstance(v, bool):
        return "on" if v else "off"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def render_report(result: ClaimResult, llm_summary: str | None = None) -> str:
    env = _env()
    tpl = env.get_template("report.html.j2")
    runs = [r.to_dict() for r in result.runs]
    key_stat_names = list(result.runs[0].key_stats.keys()) if result.runs else []
    provenance = round_floats(result.provenance)
    return tpl.render(
        spec=result.spec,
        spec_yaml=result.spec.to_yaml().strip(),
        mode=result.mode,
        verdict=result.verdict,
        verdict_icon=ICONS[result.verdict.label],
        data=result.data_summary,
        thresholds=result.thresholds,
        runs=runs,
        key_stat_names=key_stat_names,
        run_strip=run_strip_svg(runs),
        by_choice_chart=survival_by_choice_svg(result.verdict.by_choice),
        comparison_html=_comparison_html(result),
        as_labelled=result.as_labelled,
        as_labelled_html=_as_labelled_html(result),
        reference_source=result.reference_source,
        methods=result.methods,
        provenance=provenance,
        provenance_json=json.dumps(provenance, indent=2, default=json_default),
        total_runtime_s=result.total_runtime_s,
        llm_summary=llm_summary,
    )


def write_report(result: ClaimResult, path: Path, llm_summary: str | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_report(result, llm_summary), encoding="utf-8")
    return path


# ------------------------------------------------------------------ sections


def _comparison_html(result: ClaimResult) -> str:
    ct = result.spec.claim_type
    base = next((r for r in result.runs if r.is_baseline), result.runs[0] if result.runs else None)
    if base is None:
        return "<p>No runs.</p>"
    if ct == "population_exists":
        m = base.match
        rows = [
            ("Matched cells (baseline run)", m.get("n_cells")),
            ("Marker AUROC", fmtnum(m.get("marker_auroc"))),
            ("Jaccard vs reference", fmtnum(m.get("jaccard"))),
            ("Samples with ≥ min cells", m.get("n_samples_present")),
            ("Largest single-sample share", fmtnum(m.get("max_sample_share"))),
        ]
        body = "".join(
            f"<tr><td>{escape(str(k))}</td><td class='num'>{escape(str(v))}</td></tr>"
            for k, v in rows
        )
        note = (
            "<p class='small'>For a population claim the replicate-aware check is whether the "
            "population is "
            "present in several samples rather than driven by one. Recovered-as-a-cluster versus "
            "recovered-and-replicated across the grid: "
            f"{_pct(result.verdict.recovered_survival)} vs "
            f"{_pct(result.verdict.survival)}.</p>"
        )
        return f"<table><tbody>{body}</tbody></table>{note}"
    if ct == "composition_shift":
        return _composition_table(base.details, base.label) + _composition_grid_table(result)
    return _de_table(base.details, "Baseline run") + _de_grid_summary(result)


def _composition_table(d: dict[str, Any], label: str) -> str:
    rows = [
        ("Samples (A vs B)", f"{d.get('n_samples_a')} vs {d.get('n_samples_b')}"),
        ("Cells (A vs B)", f"{d.get('n_cells_a')} vs {d.get('n_cells_b')}"),
        (
            "Mean proportion A / B",
            f"{fmtnum(d.get('mean_proportion_a'))} / {fmtnum(d.get('mean_proportion_b'))}",
        ),
        ("Effect (log-odds A − B)", fmtnum(d.get("effect_log_odds_a_minus_b"))),
        ("Replicate-aware: Welch t on logit proportions, p", fmtnum(d.get("p_logit_ttest"))),
        ("Replicate-aware: Mann-Whitney U on proportions, p", fmtnum(d.get("p_mannwhitney"))),
        (
            "Cell-level: chi-square on pooled cells, p (not used)",
            fmtnum(d.get("p_cell_level_chi2")),
        ),
        ("Holds (replicate-aware)", "yes" if d.get("holds") else "no"),
        ("Holds (cell-level)", "yes" if d.get("holds_cell_level") else "no"),
    ]
    body = "".join(
        f"<tr><td>{escape(k)}</td><td class='num'>{escape(str(v))}</td></tr>" for k, v in rows
    )
    per = d.get("per_sample") or []
    per_rows = "".join(
        f"<tr><td class='mono'>{escape(str(p['sample']))}</td>"
        f"<td class='mono'>{escape(str(p['condition']))}</td>"
        f"<td class='num'>{p['k']}</td><td class='num'>{p['n']}</td>"
        f"<td class='num'>{fmtnum(p['proportion'])}</td></tr>"
        for p in per
    )
    per_tbl = (
        "<details><summary>Per-sample proportions</summary><table><thead><tr><th>Sample</th>"
        "<th>Condition</th><th class='num'>In population</th><th class='num'>Cells</th>"
        "<th class='num'>Proportion</th></tr></thead>"
        f"<tbody>{per_rows}</tbody></table></details>"
        if per
        else ""
    )
    return f"<h3>{escape(label)}</h3><table><tbody>{body}</tbody></table>{per_tbl}"


def _composition_grid_table(result: ClaimResult) -> str:
    v = result.verdict
    return (
        "<p class='small'>Across the grid the replicate-aware test held in "
        f"{_pct(v.survival)} of runs; the cell-level chi-square in "
        f"{_pct(v.cell_level_survival)}.</p>"
    )


def _de_table(d: dict[str, Any], label: str) -> str:
    per = d.get("per_gene") or []
    head = (
        "<tr><th>Gene</th><th class='num'>Cell log2FC</th><th class='num'>Cell padj</th>"
        "<th class='num'>Pseudobulk log2FC</th><th class='num'>Pseudobulk padj</th>"
        "<th class='num'>Genome-wide padj</th><th>Agreement</th></tr>"
    )
    rows = "".join(
        f"<tr><td class='mono'>{escape(str(g['gene']))}</td>"
        f"<td class='num'>{fmtnum(g.get('cell_log2fc'))}</td>"
        f"<td class='num'>{fmtnum(g.get('cell_padj'))}</td>"
        f"<td class='num'>{fmtnum(g.get('pb_log2fc'))}</td>"
        f"<td class='num'>{fmtnum(g.get('pb_padj'))}</td>"
        f"<td class='num'>{fmtnum(g.get('pb_padj_genomewide'))}</td>"
        f"<td>{_agreement(g.get('agreement'))}</td></tr>"
        for g in per
    )
    meta = (
        f"<p class='small'>{d.get('n_cells_a')} vs {d.get('n_cells_b')} cells; "
        f"{d.get('n_samples_a')} vs {d.get('n_samples_b')} pseudobulk replicates; "
        f"{d.get('n_genes_tested_pseudobulk') or '–'} genes in the DESeq2 model."
        + (f" ⚠ {'; '.join(d['warnings'])}" if d.get("warnings") else "")
        + "</p>"
    )
    return (
        f"<h3>{escape(label)}</h3><table><thead>{head}</thead><tbody>{rows}</tbody></table>{meta}"
    )


def _de_grid_summary(result: ClaimResult) -> str:
    v = result.verdict
    counts: dict[str, int] = {}
    for r in result.runs:
        for g in r.details.get("per_gene", []):
            counts[g["agreement"]] = counts.get(g["agreement"], 0) + 1
    parts = ", ".join(f"{k.replace('_', ' ')}: {n}" for k, n in sorted(counts.items()))
    return (
        f"<p class='small'>Across the grid the pseudobulk test held in {_pct(v.survival)} of runs "
        f"and the cell-level test in {_pct(v.cell_level_survival)}. "
        f"Gene × run agreement counts: {parts}.</p>"
    )


def _as_labelled_html(result: ClaimResult) -> str:
    d = result.as_labelled
    if not d:
        return ""
    ct = result.spec.claim_type
    if ct == "composition_shift":
        return _composition_table(d, "Using the labels as given")
    if ct == "differential_expression":
        return _de_table(d, "Using the labels as given")
    rows = [
        ("Labelled cells", d.get("n_cells")),
        ("Samples with ≥ min cells", d.get("n_samples_present")),
        ("Largest single-sample share", fmtnum(d.get("max_sample_share"))),
    ]
    body = "".join(
        f"<tr><td>{escape(str(k))}</td><td class='num'>{escape(str(v))}</td></tr>" for k, v in rows
    )
    return f"<table><tbody>{body}</tbody></table>"


def _agreement(a: str | None) -> str:
    if a == "both":
        return "<span class='pass'>both</span>"
    if a == "cell_only":
        return "<span class='fail'>cell-level only (pseudoreplication)</span>"
    if a == "pseudobulk_only":
        return "pseudobulk only"
    return "<span class='na'>neither</span>"


def _pct(x: float | None) -> str:
    return "–" if x is None else f"{100 * x:.0f}%"
