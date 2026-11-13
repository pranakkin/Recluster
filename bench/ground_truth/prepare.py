"""Build the single-cell inputs, compute the bulk truth and draw the claims (PLAN.md).

Inputs (downloaded into data/hagai2018/, see README.md here): the per-animal single-cell UMI
tables from E-MTAB-6754, the bulk matrices from the Squair et al. Zenodo record, and Ensembl
BioMart (release 116) gene ID -> symbol tables. Writes data/hagai2018/<species>_<stim>.h5ad,
bulk_truth/<comparison>.csv and claims/*.yaml plus claims/claims.json.

Usage: uv run python bench/ground_truth/prepare.py
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

HERE = Path(__file__).parent
DATA = Path("data/hagai2018")
COMPARISONS = [("mouse", "lps4"), ("mouse", "pic4"), ("rat", "lps4"), ("rat", "pic4")]
BIOMART = {"mouse": "mmusculus.tsv", "rat": "rnorvegicus.tsv"}
N_PER_CLASS = 10
SEED = 0


def read_sc_table(path: Path) -> tuple[list[str], list[str], sp.csr_matrix]:
    """One authors' table: first line barcodes, then one gene per line (ID then counts)."""
    with gzip.open(path, "rt") as fh:
        barcodes = fh.readline().split()
        ids, rows = [], []
        for line in fh:
            gene, _, rest = line.partition(" ")
            ids.append(gene)
            rows.append(np.array(rest.split(), dtype=np.int32))
    x = sp.csr_matrix(np.vstack(rows).T)  # cells x genes
    return barcodes, ids, x


def symbol_map(species: str) -> dict[str, str]:
    """Ensembl ID -> symbol, keeping only symbols that exactly one ID maps to."""
    m = pd.read_csv(DATA / "biomart" / BIOMART[species], sep="\t").dropna()
    m = m.drop_duplicates("Gene stable ID")
    unique = m["Gene name"].map(m["Gene name"].value_counts()) == 1
    return dict(zip(m.loc[unique, "Gene stable ID"], m.loc[unique, "Gene name"], strict=True))


def build_sc(species: str, stim: str) -> ad.AnnData:
    parts = []
    for animal in (1, 2, 3):
        for cond in (stim, "unst"):
            path = DATA / "sc" / f"{species}{animal}_{cond}_filtered_by_cell_cluster0.txt.gz"
            barcodes, ids, x = read_sc_table(path)
            obs = pd.DataFrame(
                {
                    "animal": f"{species}{animal}",
                    "condition": "stim" if cond == stim else "unst",
                    "cell_type": "phagocyte",
                },
                index=[f"{species}{animal}_{cond}_{b}" for b in barcodes],
            )
            parts.append(ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=ids)))
    ids0 = list(parts[0].var_names)
    assert all(list(p.var_names) == ids0 for p in parts), "gene order differs between files"
    adata = ad.concat(parts, join="inner")
    adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
    names = symbol_map(species)
    adata.var["ensembl_id"] = adata.var_names
    adata.var_names = [names.get(i, i) for i in adata.var_names]
    adata.var["mapped"] = [i in names for i in adata.var["ensembl_id"]]
    for c in ("animal", "condition", "cell_type"):
        adata.obs[c] = adata.obs[c].astype("category")
    return adata


def bulk_truth(species: str, stim: str) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    key = f"Hagai2018_{species}-{'lps' if stim == 'lps4' else 'pic'}"
    mdir = DATA / "bulk" / "bulk_rnaseq" / "matrices"
    counts = pd.read_csv(mdir / f"{key}_matrix.tsv.gz", sep="\t")
    counts.index = pd.read_csv(mdir / f"{key}_features.tsv.gz", sep="\t")["x"].to_numpy()
    meta = pd.read_csv(mdir / f"{key}_meta.tsv.gz", sep="\t").set_index("sample")
    meta["condition"] = np.where(meta["label"] == "UNST", "unst", "stim")
    counts = counts.loc[~counts.index.duplicated(keep=False)]
    total = counts.sum(axis=1)
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts.T.loc[meta.index].astype(int),
        metadata=meta[["condition"]],
        design="~condition",
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "stim", "unst"], inference=inference, quiet=True)
    st.summary()
    res = st.results_df[["baseMean", "log2FoldChange", "pvalue", "padj"]].copy()
    res["bulk_total"] = total.reindex(res.index)
    lfc, padj = res["log2FoldChange"], res["padj"]
    res["truth"] = np.select(
        [(padj < 0.01) & (lfc.abs() >= 1), (padj >= 0.5) & (lfc.abs() < 0.5)],
        ["change", "no_change"],
        default="",
    )
    return res


def draw(adata: ad.AnnData, truth: pd.DataFrame) -> dict[str, list[str]]:
    detected = np.asarray((adata.X > 0).mean(axis=0)).ravel()
    sc = pd.DataFrame({"detected": detected, "mapped": adata.var["mapped"].to_numpy()},
                      index=adata.var_names)  # fmt: skip
    sc = sc[sc["mapped"] & ~sc.index.duplicated(keep=False)]
    eligible = sc.index[sc["detected"] >= 0.10].intersection(truth.index[truth["bulk_total"] >= 10])
    rng = np.random.default_rng(SEED)
    out = {}
    for cls in ("change", "no_change"):
        pool = sorted(g for g in eligible if truth.at[g, "truth"] == cls)
        k = min(N_PER_CLASS, len(pool))
        out[cls] = sorted(rng.choice(pool, size=k, replace=False).tolist())
        out[f"n_{cls}_eligible"] = len(pool)  # type: ignore[assignment]
    return out


def main() -> None:
    (HERE / "bulk_truth").mkdir(exist_ok=True)
    (HERE / "claims").mkdir(exist_ok=True)
    claims = []
    for species, stim in COMPARISONS:
        comp = f"{species}_{stim}"
        h5 = DATA / f"{comp}.h5ad"
        adata = build_sc(species, stim)
        adata.write_h5ad(h5)
        truth = bulk_truth(species, stim)
        truth.to_csv(HERE / "bulk_truth" / f"{comp}.csv")
        drawn = draw(adata, truth)
        print(comp, adata.shape, {k: v if isinstance(v, int) else len(v) for k, v in drawn.items()})
        for cls in ("change", "no_change"):
            for gene in drawn[cls]:
                lfc = float(truth.at[gene, "log2FoldChange"])
                direction = "up_in_a" if lfc > 0 else "down_in_a"
                if cls == "no_change":
                    direction = "either"
                spec = {
                    "claim_type": "differential_expression",
                    "name": f"{comp} {gene} ({cls})",
                    "columns": {"condition": "condition", "sample": ["animal", "condition"],
                                "labels": "cell_type"},
                    "target": {"label_value": "phagocyte"},
                    "comparison": {"group_a": "stim", "group_b": "unst", "direction": direction},
                    "genes": [gene],
                }  # fmt: skip
                file = f"claims/{comp}_{gene}.yaml"
                (HERE / file).write_text(yaml.safe_dump(spec, sort_keys=False))
                claims.append(
                    {
                        "comparison": comp,
                        "gene": gene,
                        "truth": cls,
                        "bulk_log2fc": lfc,
                        "bulk_padj": float(truth.at[gene, "padj"]),
                        "data": str(h5),
                        "file": file,
                        "n_eligible": drawn[f"n_{cls}_eligible"],
                    }
                )
    (HERE / "claims" / "claims.json").write_text(json.dumps(claims, indent=1))
    print(len(claims), "claims")


if __name__ == "__main__":
    main()
