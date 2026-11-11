"""Inline SVG charts. No JavaScript, no external assets, verdicts never by colour alone."""

from __future__ import annotations

from html import escape
from typing import Any

PASS = "#2f6b3a"
FAIL = "#9b3b2e"
INK = "#1b1b1b"
MUTED = "#6b6b6b"
RULE = "#d8d6cf"
FONT = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"


def run_strip_svg(runs: list[dict[str, Any]]) -> str:
    """One square per run, in grid order. Filled + tick = held; hatched + cross = failed."""
    size, gap, per_row = 22, 5, 24
    n = len(runs)
    rows = (n + per_row - 1) // per_row
    w = per_row * (size + gap)
    h = rows * (size + gap)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        f'role="img" aria-label="survival strip">',
        _hatch_def(),
    ]
    for i, r in enumerate(runs):
        x = (i % per_row) * (size + gap)
        y = (i // per_row) * (size + gap)
        held = bool(r["holds"])
        fill = PASS if held else "url(#hatch)"
        stroke = INK if r["is_baseline"] else (PASS if held else FAIL)
        sw = 2.5 if r["is_baseline"] else 1
        title = escape(f"run {r['index']}: {r['label']} — {'held' if held else 'failed'}")
        out.append(
            f"<g><title>{title}</title>"
            f'<rect x="{x}" y="{y}" width="{size}" height="{size}" rx="2" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"/>'
        )
        glyph = "✓" if held else "✕"
        color = "#ffffff" if held else FAIL
        out.append(
            f'<text x="{x + size / 2}" y="{y + size / 2 + 5}" text-anchor="middle" font-size="14" '
            f'font-weight="700" fill="{color}" font-family="{FONT}">{glyph}</text></g>'
        )
    out.append("</svg>")
    return "".join(out)


def survival_by_choice_svg(by_choice: list[dict[str, Any]]) -> str:
    """For every axis value: the fraction of runs with that value in which the claim held."""
    label_w, bar_w, right_w = 250, 420, 90
    row_h, head_h, pad = 20, 26, 8
    axes: list[str] = []
    for row in by_choice:
        if row["axis"] not in axes:
            axes.append(row["axis"])
    n_rows = len(by_choice)
    h = pad * 2 + len(axes) * head_h + n_rows * row_h
    w = label_w + bar_w + right_w
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        f'role="img" aria-label="survival by analysis choice">',
        _hatch_def(),
    ]
    y = pad
    for axis in axes:
        out.append(
            f'<text x="0" y="{y + 16}" font-size="12.5" font-weight="600" fill="{MUTED}" '
            f'font-family="-apple-system, Helvetica, Arial, sans-serif" letter-spacing="0.04em">'
            f"{escape(axis.upper())}</text>"
        )
        y += head_h
        for row in (r for r in by_choice if r["axis"] == axis):
            frac = row["n_holds"] / row["n_runs"] if row["n_runs"] else 0.0
            val = _fmt(row["value"]) + ("  ●" if row["is_baseline"] else "")
            out.append(
                f'<text x="16" y="{y + 14}" font-size="13" fill="{INK}" '
                f'font-family="{FONT}">{escape(val)}</text>'
            )
            out.append(
                f'<rect x="{label_w}" y="{y + 3}" width="{bar_w}" height="{row_h - 8}" '
                'fill="url(#hatch)" '
                f'stroke="{RULE}" stroke-width="0.5"/>'
            )
            if frac > 0:
                out.append(
                    f'<rect x="{label_w}" y="{y + 3}" width="{bar_w * frac:.1f}" '
                    f'height="{row_h - 8}" fill="{PASS}"/>'
                )
            out.append(
                f'<text x="{label_w + bar_w + 10}" y="{y + 14}" font-size="13" fill="{INK}" '
                f'font-family="{FONT}">{row["n_holds"]}/{row["n_runs"]}</text>'
            )
            y += row_h
    out.append("</svg>")
    return "".join(out)


def _hatch_def() -> str:
    return (
        '<defs><pattern id="hatch" patternUnits="userSpaceOnUse" width="6" height="6" '
        'patternTransform="rotate(45)">'
        '<rect width="6" height="6" fill="#f6f4ee"/>'
        f'<line x1="0" y1="0" x2="0" y2="6" stroke="{FAIL}" '
        'stroke-width="1.2" opacity="0.7"/></pattern></defs>'
    )


def _fmt(v: Any) -> str:
    if isinstance(v, bool):
        return "on" if v else "off"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)
