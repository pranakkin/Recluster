"""Render a ClaimResult into one self-contained HTML file."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

from recluster.provenance import json_default, round_floats
from recluster.report.charts import run_strip_svg, survival_by_choice_svg
from recluster.results import ClaimResult, RunEvaluation

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
        grid=result.grid,
        runs=runs,
        key_stat_names=key_stat_names,
        run_strip=run_strip_svg(runs),
        by_choice_chart=survival_by_choice_svg(result.verdict.by_choice),
        comparison_html=_comparison_html(result),
        time_html=_time_html(result),
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
            ("Statistic that judged the candidate clusters", _basis_text(m)),
            ("Jaccard vs reference", fmtnum(m.get("jaccard"))),
            ("Samples with ≥ min cells", m.get("n_samples_present")),
            ("Largest single-sample share", fmtnum(m.get("max_sample_share"))),
            (
                "Median Scrublet score, population vs all cells (informational)",
                f"{fmtnum(m.get('doublet_score_median'))} vs "
                f"{fmtnum(m.get('doublet_score_median_all'))}",
            ),
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
        return (
            f"<table><tbody>{body}</tbody></table>{note}{_doublet_cluster_html(m)}"
            f"{_doublet_population_html(result, base)}"
        )
    if ct == "composition_shift":
        return _composition_table(base.details, base.label) + _composition_grid_table(result)
    return _de_table(base.details, "Baseline run") + _de_grid_summary(result)


def _basis_text(m: dict[str, Any]) -> str:
    """Which AUROC judged the marker candidates: stratified within-sample, overall, or mixed."""
    basis = m.get("marker_basis")
    if basis is None:
        return "labels only"
    names = {
        "within-sample": "stratified within-sample AUROC",
        "overall": "whole-dataset AUROC",
        "mixed": "mixed (see per-cluster list)",
    }
    text = names.get(str(basis), str(basis))
    judged = m.get("judged_by") or {}
    if judged and basis == "mixed":
        text += ": " + ", ".join(f"cluster {c} {b}" for c, b in judged.items())
    return text


def _composition_table(d: dict[str, Any], label: str) -> str:
    # the evaluator pairs on the donor when it can; the p-values then come from paired tests
    paired = d.get("design") == "paired"
    t_label = (
        "Replicate-aware: paired t on per-donor logit differences, p"
        if paired
        else "Replicate-aware: Welch t on logit proportions, p"
    )
    rank_label = (
        "Replicate-aware: Wilcoxon signed-rank on per-donor differences, p"
        if paired
        else "Replicate-aware: Mann-Whitney U on proportions, p"
    )
    rows = [
        ("Samples (A vs B)", f"{d.get('n_samples_a')} vs {d.get('n_samples_b')}"),
        ("Cells (A vs B)", f"{d.get('n_cells_a')} vs {d.get('n_cells_b')}"),
        (
            "Mean proportion A / B",
            f"{fmtnum(d.get('mean_proportion_a'))} / {fmtnum(d.get('mean_proportion_b'))}",
        ),
        ("Effect (log-odds A − B)", fmtnum(d.get("effect_log_odds_a_minus_b"))),
        (t_label, fmtnum(d.get("p_logit_ttest"))),
        (rank_label, fmtnum(d.get("p_mannwhitney"))),
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


def _doublet_cluster_html(m: dict[str, Any]) -> str:
    """Cluster-level doublet diagnostic of the baseline run's population (informational).

    The mixture fit asks whether the population's centroid lies between two other clusters and
    whether its library size is larger than theirs. Intermediate populations look the same, so
    the numbers are shown but no verdict depends on them.
    """
    dc = m.get("doublet_cluster")
    if not dc:
        return ""
    pair = dc.get("best_pair")
    if pair:
        w = dc.get("pair_weight")
        mix = (
            f"{pair[0]} ({fmtnum(w)}) + {pair[1]} ({fmtnum(1 - w)})"
            if w is not None
            else " + ".join(str(c) for c in pair)
        )
        rows = [
            ("Best explaining pair (weights)", mix),
            ("Residual ratio, mixture vs best single cluster", fmtnum(dc.get("residual_ratio"))),
            ("Library-size ratio, population vs mixture", fmtnum(dc.get("libsize_ratio"))),
            ("Within both cutoffs", "yes" if dc.get("flagged") else "no"),
        ]
    else:
        rows = [("Not computed", dc.get("note") or "–")]
    note = (
        "<p class='small'>Informational: no cutoff on the residual ratio or the library-size "
        "ratio separates doublet clusters from intermediate populations with margin, so this "
        "check never changes the verdict.</p>"
    )
    return "<h3>Cluster-level doublet check (baseline run)</h3>" + _kv_table(rows) + note


def _doublet_population_html(result: ClaimResult, base: RunEvaluation) -> str:
    """Clusters removed as doublet populations: the baseline run in detail, the grid in one line.

    Reads the pipeline's record in ``effective_params["doublet_clusters"]`` (present only in runs
    where the step ran) and the matched clusters' statistics the engine put in the match.
    """
    title = "<h3>Doublet-population removal (baseline run)</h3>"
    info = base.effective_params.get("doublet_clusters")
    if not info:
        body = (
            "<p class='small'>Not run in the baseline run: it runs only with doublet removal on "
            "and <span class='mono'>doublet_cluster_removal</span> enabled.</p>"
        )
        return title + body + _doublet_population_grid_line(result)
    auroc_max = info.get("auroc_max")
    stats = info.get("statistics") or {}
    removed = [str(c) for c in info.get("removed") or []]
    pop = (base.match.get("doublet_population") or {}).get("statistics") or {}
    matched = [str(c) for c in base.match.get("clusters") or []]
    if not matched:
        own = "– (population not found)"
    else:
        own = ", ".join(f"cluster {c}: {_population_stat(c, pop)}" for c in matched)
    rows: list[tuple[str, Any]] = [
        ("Matched population, AUROC vs simulated doublets", own),
        ("Removal rule", f"AUROC below {fmtval(auroc_max)}"),
        (
            "Clusters removed (cells)",
            f"{len(removed)} ({info.get('n_cells', 0)})" if removed else "none",
        ),
    ]
    out = title + _kv_table(rows)
    if removed:
        parents = info.get("parents") or {}
        body = "".join(
            f"<tr><td class='mono'>{escape(c)}</td>"
            f"<td class='num'>{fmtnum(stats.get(c))}</td>"
            f"<td class='mono'>{escape(_parent_text(parents.get(c)))}</td></tr>"
            for c in removed
        )
        out += (
            "<table><thead><tr><th>Removed cluster</th><th class='num'>AUROC</th>"
            "<th>Parent groups (clusters)</th></tr></thead>"
            f"<tbody>{body}</tbody></table>"
        )
    out += (
        "<p class='small'>Each cluster of at least 10 cells is compared with simulated doublets of "
        "two other cell-type groups; a cross-validated classifier that cannot tell them apart "
        "(low AUROC) marks a population made of doublets. Removed clusters are dropped before "
        "the population is matched, and their cells are counted with the removed doublets.</p>"
    )
    return out + _doublet_population_grid_line(result)


def _population_stat(cluster: str, stats: dict[str, Any]) -> str:
    if cluster not in stats:
        return "not tested (under 10 cells)"
    if stats[cluster] is None:
        return "not computed (fewer than two other cell-type groups)"
    return fmtnum(stats[cluster])


def _parent_text(pair: Any) -> str:
    """'1,4 + 0' for parent groups given as two lists of member clusters."""
    if not pair:
        return "–"
    return " + ".join(",".join(str(c) for c in side) for side in pair)


def _doublet_population_grid_line(result: ClaimResult) -> str:
    ran = [r for r in result.runs if r.effective_params.get("doublet_clusters")]
    if not ran:
        return ""
    hits = [
        (r.index, [str(c) for c in r.effective_params["doublet_clusters"].get("removed") or []])
        for r in ran
    ]
    hits = [(i, cl) for i, cl in hits if cl]
    detail = "; ".join(
        f"run {i}: cluster{'s' if len(cl) > 1 else ''} {', '.join(cl)}" for i, cl in hits
    )
    return (
        f"<p class='small'>Across the grid, clusters were removed as doublet populations in "
        f"{len(hits)} of the {len(ran)} runs where the step ran"
        + (f" ({escape(detail)})" if hits else "")
        + ".</p>"
    )


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


TIME_BADGE = {
    "NO TIME CONFOUND DETECTED": ("robust", "✓"),
    "POSSIBLE TIME CONFOUND": ("fragile", "◐"),
    "NOT ASSESSABLE": ("neutral", "–"),
}


def _time_html(result: ClaimResult) -> str:
    tc = result.time_check
    if not tc:
        return (
            "<p class='small'>Not applicable: this claim has no condition comparison, so there is "
            "no group contrast that collection time could confound.</p>"
        )
    label = str(tc.get("label") or "NOT ASSESSABLE")
    cls, icon = TIME_BADGE.get(label, ("neutral", "–"))
    mode = {"A": "Mode A, collection time known", "B": "Mode B, collection time unknown"}.get(
        str(tc.get("mode")), "not applicable"
    )
    reasons = "".join(f"<li>{escape(str(r))}</li>" for r in tc.get("reasons") or [])
    parts = [
        f"<div class='verdict'><div class='badge {cls} small-badge'>"
        f"<span class='icon'>{icon}</span>{escape(label)}</div>"
        f"<div class='verdict-text'><p><strong>{escape(mode)}.</strong></p>"
        f"<ul class='small'>{reasons}</ul>"
        f"<p class='small'>{escape(str(tc.get('disclaimer') or ''))}</p></div></div>"
    ]
    v = result.verdict
    # the label before the check is always recorded; only a change is worth a flag
    if v.label_before_time_check and v.label_before_time_check != v.label:
        parts.append(
            "<div class='note flag'>Verdict changed from "
            f"<strong>{escape(v.label_before_time_check)}"
            f"</strong> to <strong>{escape(v.label)}</strong> by this check.</div>"
        )
    elif v.label_before_time_check:
        parts.append("<p class='small'>This check did not change the verdict.</p>")
    bal = tc.get("balance")
    if bal and bal.get("assessable"):
        parts.append(
            "<h3>Balance of collection time across conditions</h3>" + _kv_table(_balance_rows(bal))
        )
    adj = tc.get("adjustment")
    if adj:
        parts.append(_adjustment_html(adj))
    panel = tc.get("panel")
    if panel:
        parts.append(_panel_html(panel))
    rh = tc.get("rhythmic")
    if rh:
        note = str(rh.get("note") or "")
        if not rh.get("available"):
            parts.append(f"<p class='small'>Day-night gene list: {escape(note)}.</p>")
        elif rh.get("flagged"):
            parts.append(
                "<p class='small'>Claimed genes with documented day-night variation "
                "(Wucher et al. 2023): "
                f"<span class='mono'>{escape(', '.join(str(g) for g in rh['flagged']))}</span>"
                + (f" ({escape(note)})" if note else "")
                + ".</p>"
            )
        else:
            parts.append(
                "<p class='small'>None of the claimed genes is on the public day-night list"
                + (f" ({escape(note)})" if note else "")
                + ".</p>"
            )
    return "".join(parts)


def _permutation_text(p: Any, n: Any, exact: Any, min_p: Any) -> str:
    """p with how it was obtained: exact enumeration or seeded draws, and the floor it can reach."""
    detail = []
    if n is not None:
        detail.append(f"{'exact, ' if exact else ''}{n} permutations")
    if min_p is not None:
        detail.append(f"smallest attainable p {fmtnum(min_p)}")
    return fmtnum(p) + (f" ({'; '.join(detail)})" if detail else "")


def _units_row(d: dict[str, Any]) -> list[tuple[str, Any]]:
    """The level the test ran at (samples, or donors when samples nest in donors)."""
    col = d.get("unit_column")
    return [("Units", f"{col} levels")] if col is not None else []


def _balance_rows(bal: dict[str, Any]) -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = _units_row(bal)
    if bal.get("n_samples_a") is not None:
        rows.append(("Samples (A vs B)", f"{bal.get('n_samples_a')} vs {bal.get('n_samples_b')}"))
    if bal.get("n_units_a") is not None:
        rows.append(("Units (A vs B)", f"{bal.get('n_units_a')} vs {bal.get('n_units_b')}"))
    if bal.get("circular_mean_hours_a") is not None:
        rows.append(
            (
                "Circular mean collection time A / B (h)",
                f"{bal.get('circular_mean_hours_a')} / {bal.get('circular_mean_hours_b')}",
            )
        )
    diff = bal.get("circular_mean_difference_hours")
    if diff is not None:
        rows.append(("Difference (h)", f"{float(diff):+.2f}"))
    if bal.get("mean_vector_distance") is not None:
        rows.append(("Mean-vector distance", fmtnum(bal.get("mean_vector_distance"))))
    rows.append(
        (
            "Permutation p (time vs condition)",
            _permutation_text(
                bal.get("p_permutation"),
                bal.get("n_permutations"),
                bal.get("exact"),
                bal.get("min_attainable_p"),
            ),
        )
    )
    if bal.get("p_cos_ttest") is not None or bal.get("p_sin_ttest") is not None:
        rows.append(
            (
                "Welch t on cos / sin, p",
                f"{fmtnum(bal.get('p_cos_ttest'))} / {fmtnum(bal.get('p_sin_ttest'))}",
            )
        )
    if bal.get("imbalanced") is not None:
        basis = bal.get("imbalance_basis")
        rows.append(
            (
                "Time imbalanced between conditions",
                _yes_no(bal.get("imbalanced")) + (f" ({basis})" if basis else ""),
            )
        )
    return rows


def _adjustment_html(adj: dict[str, Any]) -> str:
    rows: list[tuple[str, Any]] = []
    if adj.get("design") is not None:
        rows.append(("Design", adj.get("design")))
    if adj.get("estimable") is False:
        reason = adj.get("reason") or adj.get("estimable_reason")
        rows.append(("Estimable", "no" + (f": {reason}" if reason else "")))
    per_gene = adj.get("per_gene") or []
    # a composition refit carries one coefficient; a DE refit carries per-gene rows
    composition = "effect_unadjusted" in adj
    if composition:
        rows += [
            (
                "Effect (log-odds) unadjusted / adjusted",
                f"{fmtnum(adj.get('effect_unadjusted'))} / {fmtnum(adj.get('effect_adjusted'))}",
            ),
            (
                "p unadjusted / adjusted",
                f"{fmtnum(adj.get('p_unadjusted'))} / {fmtnum(adj.get('p_adjusted'))}",
            ),
            ("Shrinkage", fmtnum(adj.get("shrinkage"))),
        ]
    if adj.get("held_unadjusted") is not None:
        rows.append(("Held before adjustment", _yes_no(adj.get("held_unadjusted"))))
    rows.append(("Survives adjustment", _survives(adj.get("survives_adjustment"))))
    if adj.get("shrunk") is not None:
        rows.append(("Effect shrunk", _yes_no(adj.get("shrunk"))))
    if adj.get("flipped") is not None:
        rows.append(("Direction flipped", _yes_no(adj.get("flipped"))))
    title = (
        "Composition refitted with time covariates"
        if composition
        else "Claim refitted with cos(2πt/24) and sin(2πt/24) as covariates"
    )
    out = f"<h3>{title}</h3>" + _kv_table(rows)
    if per_gene:
        head = (
            "<tr><th>Gene</th><th class='num'>log2FC unadjusted</th><th class='num'>padj</th>"
            "<th class='num'>log2FC adjusted for time</th><th class='num'>padj</th>"
            "<th class='num'>Shrinkage</th><th>Survives</th></tr>"
        )
        body = "".join(
            f"<tr><td class='mono'>{escape(str(g.get('gene')))}</td>"
            f"<td class='num'>{fmtnum(g.get('lfc_unadjusted'))}</td>"
            f"<td class='num'>{fmtnum(g.get('padj_unadjusted'))}</td>"
            f"<td class='num'>{fmtnum(g.get('lfc_adjusted'))}</td>"
            f"<td class='num'>{fmtnum(g.get('padj_adjusted'))}</td>"
            f"<td class='num'>{fmtnum(g.get('shrinkage'))}</td>"
            f"<td>{_survives_span(g.get('survives'))}"
            f"{' (direction flips)' if g.get('direction_flipped') else ''}</td></tr>"
            for g in per_gene
        )
        out += f"<table><thead>{head}</thead><tbody>{body}</tbody></table>"
    return out


def _panel_html(panel: dict[str, Any]) -> str:
    rows: list[tuple[str, Any]] = [
        ("Panel", ", ".join(str(g) for g in panel.get("panel") or [])),
        ("Source", panel.get("panel_citation") or "–"),
        (
            "Genes present / usable",
            f"{len(panel.get('genes_present') or [])} / {len(panel.get('genes_usable') or [])} "
            f"(≥ {panel.get('min_total_counts')} counts)",
        ),
    ]
    if panel.get("assessable"):
        rows += _units_row(panel)
        if panel.get("n_samples_a") is not None:
            rows.append(
                ("Samples (A vs B)", f"{panel.get('n_samples_a')} vs {panel.get('n_samples_b')}")
            )
        if panel.get("index_genes"):
            rows.append(("Genes in the index", ", ".join(str(g) for g in panel["index_genes"])))
        rows += [
            (
                "Clock-state index difference (A − B)",
                fmtnum(panel.get("index_difference_a_minus_b")),
            ),
            ("Cohen's d", fmtnum(panel.get("index_cohens_d"))),
            (
                "Permutation p on the index",
                _permutation_text(
                    panel.get("index_p_permutation"),
                    panel.get("index_n_permutations"),
                    panel.get("index_exact"),
                    panel.get("min_attainable_p"),
                ),
            ),
        ]
        if panel.get("index_welch_p") is not None:
            rows.append(("Welch t on the index, p", fmtnum(panel.get("index_welch_p"))))
        rows.append(("Panel genes differing (BH < α)", str(panel.get("n_panel_genes_different"))))
    else:
        rows.append(("Not assessable", panel.get("reason") or ""))
    out = "<h3>Clock-gene state across samples</h3>" + _kv_table(rows)
    per_gene = panel.get("per_gene") or []
    if panel.get("assessable") and per_gene:
        head = (
            "<tr><th>Gene</th><th>Arm</th><th class='num'>mean log2CPM A</th>"
            "<th class='num'>mean log2CPM B</th><th class='num'>A − B</th>"
            "<th class='num'>padj</th></tr>"
        )
        body = "".join(
            f"<tr><td class='mono'>{escape(str(g.get('gene')))}</td>"
            f"<td>{escape(str(g.get('arm') or ''))}</td>"
            f"<td class='num'>{fmtnum(g.get('mean_log2cpm_a'))}</td>"
            f"<td class='num'>{fmtnum(g.get('mean_log2cpm_b'))}</td>"
            f"<td class='num'>{fmtnum(g.get('diff_a_minus_b'))}</td>"
            f"<td class='num'>{fmtnum(g.get('padj'))}</td></tr>"
            for g in per_gene
        )
        out += (
            "<details><summary>Per-gene panel values</summary><table>"
            f"<thead>{head}</thead><tbody>{body}</tbody></table></details>"
        )
    return out


def _yes_no(v: Any) -> str:
    return "–" if v is None else ("yes" if v else "no")


def _survives(v: Any) -> str:
    """Three-valued: the adjusted model may not be estimable (None)."""
    return "not estimable" if v is None else ("yes" if v else "no")


def _survives_span(v: Any) -> str:
    if v is None:
        return "<span class='na'>not estimable</span>"
    return f"<span class='{'pass' if v else 'fail'}'>{'yes' if v else 'no'}</span>"


def _kv_table(rows: list[tuple[str, Any]]) -> str:
    body = "".join(
        f"<tr><td>{escape(str(k))}</td><td class='num'>{escape(str(v))}</td></tr>" for k, v in rows
    )
    return f"<table><tbody>{body}</tbody></table>"
