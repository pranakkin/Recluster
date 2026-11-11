"""Build recluster/resources/daynight_genes.csv from Wucher et al. 2023 (PLoS Biol, CC BY 4.0).

Source: S8 Table of Wucher V, Sousa M, Bonnet A, Vauclair A, Dupuy A, Dabout A, et al.
Day-night and seasonal variation of human gene expression across tissues. PLoS Biol
2023;21(2):e3001986. https://doi.org/10.1371/journal.pbio.3001986.s030

The table lists, for every day-night gene, the number of GTEx tissues in which it is higher
during the day or during the night (all tissues, non-brain, brain) with a binomial P value.
It is gene-level: per-tissue membership is not in the public tables, so Recluster's overlap
check is tissue-agnostic and says so.

Usage: python scripts/fetch_daynight_genes.py   (needs openpyxl, in the dev group)
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

import pandas as pd

URL = "https://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.3001986.s030&type=supplementary"
OUT = Path(__file__).resolve().parents[1] / "recluster" / "resources" / "daynight_genes.csv"


def main() -> None:
    cache = Path("data") / "wucher2023_S8.xlsx"
    cache.parent.mkdir(exist_ok=True)
    if not cache.exists():
        req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r:  # noqa: S310 - fixed https URL
            cache.write_bytes(r.read())
    print(f"S8 Table sha256 {hashlib.sha256(cache.read_bytes()).hexdigest()[:16]}…")
    df = pd.read_excel(cache, sheet_name=0)
    out = pd.DataFrame(
        {
            "gene": df["Name"].astype(str),
            "ensembl": df["Gene"].astype(str).str.split(".").str[0],
            "day_tissues": df["Day"].astype(int),
            "night_tissues": df["Night"].astype(int),
            "total_tissues": df["Total"].astype(int),
            "p_binomial_all_tissues": df["pvalue all tissues"],
        }
    )
    out = out[out["total_tissues"] >= 1].sort_values(
        ["total_tissues", "gene"], ascending=[False, True]
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Day-night genes from Wucher et al. 2023, PLoS Biol 21(2):e3001986, S8 Table "
        "(CC BY 4.0).\n"
        "# Columns: gene symbol, Ensembl id, number of GTEx tissues where the gene is higher "
        "by day,\n"
        "# by night, total, binomial P over all tissues. Gene-level: no per-tissue membership.\n"
    )
    OUT.write_text(header + out.to_csv(index=False))
    print(f"wrote {OUT}: {len(out)} genes; >=2 tissues: {(out['total_tissues'] >= 2).sum()}")


if __name__ == "__main__":
    main()
