"""Plain-English summary of results the code already computed.

The model sees a short digest of numbers and is told to restate them, nothing more. A guard
rejects any summary that contains a number absent from the digest, so the text can never
introduce a figure the code did not produce. The output is always labelled as LLM output.
"""

from __future__ import annotations

from typing import Any

from recluster.llm.common import get_client, model_name, numbers_in, text_of

LABEL = "LLM summary of computed results"

SYSTEM = (
    "You write short plain-English summaries of a computational result for a biologist. "
    "You are given a digest of numbers computed by deterministic code. Rules: use only numbers "
    "that appear in the digest, copied exactly; do not compute new numbers, percentages or "
    "differences; do not speculate about biology or mechanism beyond what the digest states; "
    "do not give recommendations. Write 3 to 5 sentences in one paragraph: the verdict and what "
    "it means, the survival across the grid, the main setting(s) that broke the claim if any, "
    "and the replicate-aware versus cell-level comparison when present. If the digest says the "
    "pseudoreplication flag is set, say in one sentence that cell-level statistics treat cells as "
    "independent replicates and overstate significance. Plain prose, no bullet points, no headings."
)


def digest_from_provenance(p: dict[str, Any]) -> str:
    """A compact, deterministic text digest; every number the summary may use is in here."""
    v = p["verdict"]
    claim = p["claim"]
    n, k = v["n_runs"], v["n_holds"]
    lines = [
        f"claim: {claim['name']}",
        f"claim_type: {claim['claim_type']}",
        f"verdict: {v['label']}",
        f"runs_held: {k} of {n} runs; runs_failed: {n - k}; "
        f"survival_percent: {round(100 * v['survival'])}",
        f"baseline_holds: {'yes' if v['baseline_holds'] else 'no'}",
        f"thresholds: ROBUST at or above {round(100 * p['thresholds']['robust_min_survival'])} "
        "percent survival; "
        f"FRAGILE at or above {round(100 * p['thresholds']['fragile_min_survival'])} percent; "
        f"significance alpha {p['thresholds']['alpha']:g}",
    ]
    if v.get("cell_level_survival") is not None and claim["claim_type"] != "population_exists":
        lines.append(
            f"cell_level_survival_percent (naive, not used for the verdict): "
            f"{round(100 * v['cell_level_survival'])}"
        )
    lines.append(
        f"pseudoreplication_flag: {'set' if v.get('pseudoreplication_flag') else 'not set'}"
    )
    if v.get("oat_flips"):
        flips = "; ".join(
            ", ".join(f"{a}={_fmt(b)}" for a, b in f["params_changed"].items())
            for f in v["oat_flips"][:6]
        )
        verb = "broke" if v["baseline_holds"] else "rescued"
        lines.append(f"single_changes_that_{verb}_it: {flips}")
    if v.get("breakers"):
        top = "; ".join(
            f"{b['axis']}={_fmt(b['value'])} failed in {b['n_runs'] - b['n_holds']} "
            f"of {b['n_runs']} runs"
            for b in v["breakers"][:4]
        )
        lines.append(f"failure_rates_by_setting: {top}")
    base = next((r for r in p["runs"] if r.get("is_baseline")), None)
    if base:
        m = base.get("match", {})
        lines.append(
            f"baseline_population: {m.get('n_cells')} cells"
            + (
                f", marker AUROC {_fmt(m.get('marker_auroc'))}"
                if m.get("marker_auroc") is not None
                else ""
            )
            + (f", Jaccard {_fmt(m.get('jaccard'))}" if m.get("jaccard") is not None else "")
            + (
                f", present in {m.get('n_samples_present')} samples"
                if m.get("n_samples_present") is not None
                else ""
            )
        )
        d = base.get("details", {})
        if claim["claim_type"] == "composition_shift":
            lines.append(
                "baseline_composition: mean proportion group_a "
                f"{_fmt(d.get('mean_proportion_a'))} vs group_b "
                f"{_fmt(d.get('mean_proportion_b'))}; "
                f"replicate-aware p {_fmt(d.get('p_logit_ttest'))}; "
                f"cell-level chi-square p {_fmt(d.get('p_cell_level_chi2'))}; "
                f"samples {d.get('n_samples_a')} vs {d.get('n_samples_b')}"
            )
        if claim["claim_type"] == "differential_expression":
            for g in d.get("per_gene", [])[:6]:
                lines.append(
                    f"baseline_gene {g['gene']}: pseudobulk log2FC {_fmt(g.get('pb_log2fc'))}, "
                    f"pseudobulk padj {_fmt(g.get('pb_padj'))}, "
                    f"cell-level padj {_fmt(g.get('cell_padj'))}, "
                    f"agreement {g.get('agreement')}"
                )
    al = p.get("as_labelled")
    if al and "holds" in al:
        lines.append(
            f"as_labelled (using the existing labels without re-clustering): "
            f"replicate-aware holds {'yes' if al['holds'] else 'no'}; "
            f"cell-level holds {'yes' if al.get('holds_cell_level') else 'no'}"
        )
    return "\n".join(lines)


def _fmt(x: Any) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "on" if x else "off"
    if isinstance(x, float):
        if x != 0 and abs(x) < 1e-3:
            return f"{x:.1e}"
        return f"{x:.3g}"
    return str(x)


class SummaryRejected(RuntimeError):
    """The model introduced a number that is not in the digest."""


def check_numbers(summary: str, digest: str) -> set[str]:
    """Numbers in ``summary`` that do not appear in ``digest`` (empty set = ok)."""
    allowed = numbers_in(digest)
    return {n for n in numbers_in(summary) if n not in allowed}


def summarize_provenance(
    provenance: dict[str, Any],
    *,
    client: Any = None,
    model: str | None = None,
    max_attempts: int = 2,
) -> tuple[str, str]:
    """Return (summary_text, model_id). Raises SummaryRejected if the guard fails twice."""
    client = client or get_client()
    model_id = model_name(model)
    digest = digest_from_provenance(provenance)
    user = (
        "Digest of computed results (the only source of numbers you may use):\n\n"
        f"{digest}\n\nWrite the summary now."
    )
    messages: list[dict[str, Any]] = [{"role": "user", "content": user}]
    for _ in range(max_attempts):
        response = client.messages.create(
            model=model_id,
            max_tokens=1024,
            system=SYSTEM,
            messages=messages,
        )
        if getattr(response, "stop_reason", "") == "refusal":
            raise SummaryRejected("the model declined to write a summary")
        text = text_of(response)
        bad = check_numbers(text, digest)
        if not bad:
            return text, model_id
        messages += [
            {"role": "assistant", "content": text},
            {
                "role": "user",
                "content": (
                    f"Your summary used numbers not in the digest: {sorted(bad)}. Rewrite it using "
                    "only numbers copied exactly from the digest."
                ),
            },
        ]
    raise SummaryRejected(f"summary kept introducing numbers not in the digest: {sorted(bad)}")


def summary_markdown(summary: str, model_id: str) -> str:
    return (
        f"## {LABEL}\n\n_Written by {model_id} from the numbers in provenance.json; "
        f"not a computation._\n\n{summary}\n"
    )
