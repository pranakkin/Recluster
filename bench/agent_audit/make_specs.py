"""Turn the agent's stated claims into claim specs, following PROTOCOL.md and its amendments.

Every quote must appear verbatim in run/final_answer.md (checked below), so the mapping from the
agent's words to a spec can be read off this file. Where a statement names a group of cell types,
the members are computed from the agent's own table (amendment 2); the remaining members of the
group are tested as well, in specs marked secondary (amendment 6). Writes claims/*.yaml and
claims/claims.json.

Usage: uv run python bench/agent_audit/make_specs.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).parent
ANSWER = (HERE / "run" / "final_answer.md").read_text()
AGENT_DE = pd.read_csv(HERE / "run" / "workspace" / "pseudobulk_DE_stim_vs_ctrl.csv")
OUT = HERE / "claims"

COLUMNS = {"condition": "label", "sample": ["replicate", "label"], "labels": "cell_type"}
STIM_VS_CTRL = {"group_a": "stim", "group_b": "ctrl"}

# the agent's table names -> the labels it wrote to annotated.h5ad (obs['cell_type'])
POPULATIONS = {
    "| CD4 T naive (2 clusters merged) |": "CD4 T naive",
    "| CD14 monocytes |": "CD14 Mono",
    "| CD4 T memory |": "CD4 T memory",
    "| B |": "B",
    "| NK |": "NK",
    "| CD8 T |": "CD8 T",
    "| CD16 monocytes |": "CD16 Mono",
    "| T cells with high stress genes |": "T (stress/HSP-high)",
    "| cDC |": "cDC",
    "| pDC |": "pDC",
    "| Megakaryocyte/platelet |": "Megakaryocyte/platelet",
}
HEDGED = {  # labels the agent itself called uncertain ("Two labels are less certain")
    "CD4 T naive": "**CD4 T cells:** the CD4 gene itself is barely detected",
    "CD4 T memory": "**CD4 T cells:** the CD4 gene itself is barely detected",
    "Megakaryocyte/platelet": "so it probably includes platelets stuck to T cells",
}
T_TYPES = ["CD4 T naive", "CD4 T memory", "CD8 T", "T (stress/HSP-high)"]
MYELOID = ["CD14 Mono", "CD16 Mono", "cDC"]  # "monocytes and DCs/cDC" among the 9 tested types
TESTED_DE = ["B", "CD14 Mono", "CD16 Mono", *T_TYPES[:2], "CD8 T", "NK", T_TYPES[3], "cDC"]

Q_POP = "I annotated 11 cell types."
Q_COMP = "Nothing passes FDR < 0.05."
Q_SHARED = (
    "70 genes go up significantly in all 9 tested types. Examples are ISG15, IFI6, "
    "IFIT1/2/3, MX1/MX2, IRF7, OAS1, LY6E, IFITM1/2/3, RSAD2, CMPK2, HERC5 and BST2."
)
SHARED = [
    "ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "IRF7", "OAS1",
    "LY6E", "IFITM1", "IFITM2", "IFITM3", "RSAD2", "CMPK2", "HERC5", "BST2",
]  # fmt: skip
Q_MYELOID_UP = (
    "Examples are CCL8 (log2 fold change 9.8 in CD14 monocytes), CXCL10 (9.7 in CD14 "
    "monocytes, 12.0 in cDC), CXCL11, TNFSF10 and IFI27 (7.5), plus IL27, IDO1 and CH25H in "
    "CD16 monocytes."
)
MYELOID_UP = ["CCL8", "CXCL10", "CXCL11", "TNFSF10", "IFI27"]
CD16_UP = ["IL27", "IDO1", "CH25H"]
Q_MYELOID_DOWN = (
    "Monocytes and cDC lose classical-monocyte and inflammatory genes: CD14 (−3.1 in CD16 "
    "monocytes), VCAN (−4.4), CLEC5A (−4.8), MMP9 (−4.3), IL8 (−4.1), CD93 and S100A8."
)
MYELOID_DOWN = ["CD14", "VCAN", "CLEC5A", "MMP9", "IL8", "CD93", "S100A8"]
Q_NK_DOWN = (
    "Lymphocytes lose only a few genes, and modestly, e.g. AMICA1, NELL2, ITM2C and S100A4 in NK."
)
NK_DOWN = ["AMICA1", "NELL2", "ITM2C", "S100A4"]
Q_CAVEAT = 'CXCL10, CCL8 and CXCL11 also appear strongly "up" in T, B and NK cells.'
CAVEAT_UP = ["CXCL10", "CCL8", "CXCL11"]
LYMPHOID = [*T_TYPES, "B", "NK"]


def agent_supports(label: str, gene: str, direction: str) -> bool:
    """True if the agent's own table has the gene at padj < 0.05 in the stated direction."""
    row = AGENT_DE[(AGENT_DE.cell_type == label) & (AGENT_DE.gene == gene)]
    if row.empty:
        return False
    r = row.iloc[0]
    sign = r.log2FoldChange > 0 if direction == "up_in_a" else r.log2FoldChange < 0
    return bool(r.padj < 0.05 and sign)


def claims() -> list[dict]:
    out: list[dict] = []
    for row, label in POPULATIONS.items():
        out.append(
            {
                "id": f"P-{label}",
                "statement": "annotation",
                "quote": [Q_POP, row] + ([HEDGED[label]] if label in HEDGED else []),
                "agent_says": "exists" + (" (hedged by the agent)" if label in HEDGED else ""),
                "secondary": False,
                "spec": {
                    "claim_type": "population_exists",
                    "name": f"agent: {label} is a distinct population",
                    "columns": {"labels": "cell_type"},
                    "target": {"label_value": label},
                },
            }
        )
    for row, label in POPULATIONS.items():  # amendment 9: the same, with the sample column
        out.append(
            {
                "id": f"P9-{label}",
                "statement": "annotation",
                "quote": [Q_POP, row] + ([HEDGED[label]] if label in HEDGED else []),
                "agent_says": "exists" + (" (hedged by the agent)" if label in HEDGED else ""),
                "secondary": False,
                "supersedes": f"P-{label}",
                "spec": {
                    "claim_type": "population_exists",
                    "name": f"agent: {label} is a distinct population (per-sample doublets)",
                    "columns": {"sample": COLUMNS["sample"], "labels": "cell_type"},
                    "target": {"label_value": label},
                },
            }
        )
    for label in POPULATIONS.values():
        out.append(
            {
                "id": f"C-{label}",
                "statement": "abundance null",
                "quote": [Q_COMP],
                "agent_says": "no significant change after FDR over the 11 types",
                "secondary": False,
                "spec": {
                    "claim_type": "composition_shift",
                    "name": f"agent null: {label} proportion differs stim vs ctrl",
                    "columns": COLUMNS,
                    "target": {"label_value": label},
                    "comparison": {**STIM_VS_CTRL, "direction": "different"},
                },
            }
        )
    for label in TESTED_DE:
        out.append(_de("S", "shared ISGs", [Q_SHARED], label, SHARED, "up_in_a"))
    out += _group("MU", "myeloid up", [Q_MYELOID_UP], MYELOID, MYELOID_UP, "up_in_a")
    out.append(_de("MU16", "CD16 up", [Q_MYELOID_UP], "CD16 Mono", CD16_UP, "up_in_a"))
    out += _group("MD", "myeloid down", [Q_MYELOID_DOWN], MYELOID, MYELOID_DOWN, "down_in_a")
    for gene in NK_DOWN:  # one spec per gene: two of them are not in the agent's NK table
        out.append(_de(f"ND-{gene}", "NK down", [Q_NK_DOWN], "NK", [gene], "down_in_a"))
    out += _group("CV", "caveat up", [Q_CAVEAT], LYMPHOID, CAVEAT_UP, "up_in_a")
    return out


def _group(
    prefix: str, what: str, quote: list[str], members: list[str], genes: list[str], direction: str
) -> list[dict]:
    """Amendment 2 (primary: where the agent's table supports it) and 6 (secondary: the rest)."""
    out = []
    for label in members:
        backed = [g for g in genes if agent_supports(label, g, direction)]
        rest = [g for g in genes if g not in backed]
        if backed:
            out.append(_de(prefix, what, quote, label, backed, direction))
        if rest:
            literal = what + " (literal reading)"
            out.append(_de(prefix + "x", literal, quote, label, rest, direction, True))
    return out


def _de(
    prefix: str,
    what: str,
    quote: list[str],
    label: str,
    genes: list[str],
    direction: str,
    secondary: bool = False,
) -> dict:
    verb = "up" if direction == "up_in_a" else "down"
    return {
        "id": f"{prefix}-{label}",
        "statement": what,
        "quote": quote,
        "agent_says": f"{', '.join(genes)} {verb} in {label}",
        "secondary": secondary,
        "spec": {
            "claim_type": "differential_expression",
            "name": f"agent: {what} in {label}" + (f" ({genes[0]})" if len(genes) == 1 else ""),
            "columns": COLUMNS,
            "target": {"label_value": label},
            "comparison": {**STIM_VS_CTRL, "direction": direction},
            "genes": genes,
            "require": "all",
        },
    }


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("*.yaml"):
        old.unlink()
    items = claims()
    ids = [c["id"] for c in items]
    assert len(ids) == len(set(ids)), "duplicate claim ids"
    for c in items:
        for q in c["quote"]:
            if q not in ANSWER:
                raise SystemExit(f"{c['id']}: quote not found verbatim in final_answer.md: {q!r}")
        safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in c["id"])
        c["file"] = f"claims/{safe}.yaml"
        text = yaml.safe_dump(c["spec"], sort_keys=False, allow_unicode=True)
        (HERE / c["file"]).write_text(text)
    (OUT / "claims.json").write_text(json.dumps(items, indent=1, ensure_ascii=False))
    n2 = sum(c["secondary"] for c in items)
    print(f"{len(items)} specs ({len(items) - n2} primary, {n2} secondary) written to {OUT}")


if __name__ == "__main__":
    main()
