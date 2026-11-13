"""Recompute the post hoc numbers in RESULTS.md (those marked †).

The post hoc analysis looks past the 40 scored claims. It recomputes:

- the bulk DESeq2 of each LPS comparison, unpaired (the truth, with its lfcSE) and paired by
  replicate index (an animal order the plan declined to assume);
- the paired and unpaired pseudobulk DESeq2 on all cells, one sample per animal and condition,
  which is Recluster's replicate-level test before doublet removal;
- per-animal log2(stim/unst) of the three rat false calls with DESeq2 size factors and with
  library-size (CP10k) scaling, and the offset between the two scalings;
- a cycling-cell flag (at least 4 G2M or 4 S-phase marker genes detected), the share of cycling
  cells per sample and each gene's enrichment in cycling unstimulated cells, as
  log2((CP10k cycling + 0.01) / (CP10k non-cycling + 0.01));
- the one-gene paired and unpaired pseudobulk calls (p < 0.05) on every eligible gene;
- library sizes, genes detected and the rat mitochondrial share.

Inputs are data/hagai2018/<species>_lps4.h5ad and the bulk matrices (see prepare.py), plus
bulk_truth/*.csv and claims/claims.json. Writes small CSV files to posthoc/ and
posthoc/numbers.json, which maps a key for each † number to its value.

The script stops if its draw, eligibility pools or unpaired bulk DESeq2 differ from prepare.py
and bulk_truth/*.csv. When the per-claim provenance.json files from run_claims.py are present
(git ignores them), it also stops if its all-cells paired pseudobulk for the three rat false
calls differs from Recluster's own; without them that check is skipped and the outputs are the
same.

Usage: uv run python bench/ground_truth/posthoc.py
"""

from __future__ import annotations

import contextlib
import io
import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import prepare
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy import stats

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
DATA = ROOT / "data" / "hagai2018"
OUT = HERE / "posthoc"
SPECIES = ("rat", "mouse")
ANIMALS = (1, 2, 3)
FALSE_CALLS = ("Alkbh7", "Ezh2", "Son")
CYCLE_GENES = ("Top2a", "Mki67", "Cdk1", "Ube2c", "Ccnb1")
G2M = ["Mki67", "Top2a", "Cdk1", "Ube2c", "Cdca3", "Ccnb1", "Ccnb2", "Birc5", "Cenpf", "Tpx2",
       "Nusap1", "Prc1", "Hmmr", "Cdc20", "Aurka", "Aurkb", "Smc4", "Kif20a"]  # fmt: skip
S_PHASE = ["Pcna", "Mcm2", "Mcm3", "Mcm5", "Mcm6", "Rrm2", "Tyms", "Hells", "Gins2", "Ung", "Dut",
           "Uhrf1", "Clspn", "Rrm1"]  # fmt: skip
MIN_MARKERS = 4
# keeps the enrichment finite for genes absent from one group
PSEUDOCOUNT = 0.01
INFERENCE = DefaultInference(n_cpus=4)

# prepare.py reads its inputs relative to the working directory
prepare.DATA = DATA


def deseq(counts: pd.DataFrame, meta: pd.DataFrame, design: str) -> tuple[pd.DataFrame, pd.Series]:
    """DESeq2 stim against unst; returns the results table and the size factors."""
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(
            counts=counts, metadata=meta, design=design, inference=INFERENCE, quiet=True
        )
        dds.deseq2()
        st = DeseqStats(
            dds, contrast=["condition", "stim", "unst"], inference=INFERENCE, quiet=True
        )
        st.summary()
    return st.results_df, pd.Series(dds.obs["size_factors"].to_numpy(), index=dds.obs_names)


def load_bulk(species: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bulk counts (samples x genes) and metadata, read as prepare.bulk_truth reads them."""
    key = f"Hagai2018_{species}-lps"
    mdir = DATA / "bulk" / "bulk_rnaseq" / "matrices"
    counts = pd.read_csv(mdir / f"{key}_matrix.tsv.gz", sep="\t")
    counts.index = pd.read_csv(mdir / f"{key}_features.tsv.gz", sep="\t")["x"].to_numpy()
    meta = pd.read_csv(mdir / f"{key}_meta.tsv.gz", sep="\t").set_index("sample")
    meta["condition"] = np.where(meta["label"] == "UNST", "unst", "stim")
    meta["pair"] = "p" + meta["replicate"].astype(str)
    counts = counts.loc[~counts.index.duplicated(keep=False)]
    return counts.T.loc[meta.index].astype(int), meta


def sample_key(adata: ad.AnnData) -> np.ndarray:
    return (adata.obs["animal"].astype(str) + "|" + adata.obs["condition"].astype(str)).to_numpy()


def pseudobulk(x: sp.csr_matrix, names: pd.Index, key: np.ndarray) -> pd.DataFrame:
    """Summed counts per animal and condition over all cells, unique names with >= 10 counts."""
    samples = sorted(set(key))
    pb = pd.DataFrame(
        np.vstack([np.asarray(x[key == s].sum(axis=0)).ravel() for s in samples]),
        index=samples,
        columns=names,
    )
    pb = pb.loc[:, ~pb.columns.duplicated(keep=False)]
    return pb.loc[:, pb.sum(axis=0) >= 10].round().astype(int)


def eligible(adata: ad.AnnData, truth: pd.DataFrame) -> pd.Index:
    """The genes prepare.draw draws from, both classes."""
    detected = np.asarray((adata.X > 0).mean(axis=0)).ravel()
    sc = pd.DataFrame({"detected": detected, "mapped": adata.var["mapped"].to_numpy()},
                      index=adata.var_names)  # fmt: skip
    sc = sc[sc["mapped"] & ~sc.index.duplicated(keep=False)]
    return sc.index[sc["detected"] >= 0.10].intersection(truth.index[truth["bulk_total"] >= 10])


def per_animal(values: pd.Series, species: str) -> list[float]:
    """log2(stim/unst) for each animal from per-sample values indexed '<animal>|<condition>'."""
    return [
        float(np.log2(values[f"{species}{a}|stim"] / values[f"{species}{a}|unst"])) for a in ANIMALS
    ]


def cp10k(x: sp.csr_matrix, lib: np.ndarray, key: np.ndarray, mask: np.ndarray) -> pd.DataFrame:
    """Per-sample summed counts per 10,000 summed UMI over the cells in mask (genes as columns)."""
    samples = sorted(set(key))
    rows = [np.asarray(x[mask & (key == s)].sum(axis=0)).ravel() / lib[mask & (key == s)].sum()
            for s in samples]  # fmt: skip
    return pd.DataFrame(1e4 * np.vstack(rows), index=samples)


def unique_present(genes: list[str], names: pd.Index) -> list[str]:
    counts = names.value_counts()
    return [g for g in genes if counts.get(g, 0) == 1]


def n_detected(x: sp.csr_matrix, names: pd.Index, genes: list[str]) -> np.ndarray:
    return np.asarray((x[:, [names.get_loc(g) for g in genes]] > 0).sum(axis=1)).ravel()


def provenance(species: str, gene: str) -> dict | None:
    """Recluster's own all-cells ('as labelled') pseudobulk for a claim, when the run is here."""
    summary = json.loads((HERE / "out" / "summary.json").read_text())
    entry = summary.get(f"{species}_lps4/{gene}")
    path = HERE / entry["folder"] / "provenance.json" if entry else None
    if path is None or not path.exists():
        return None
    return json.loads(path.read_text())["as_labelled"]["per_gene"][0]


def analyse(species: str, claims: list[dict], num: dict, tables: dict) -> None:
    truth = pd.read_csv(HERE / "bulk_truth" / f"{species}_lps4.csv", index_col=0)
    adata = ad.read_h5ad(DATA / f"{species}_lps4.h5ad")
    x = sp.csr_matrix(adata.X)
    names = pd.Index(adata.var_names)
    key = sample_key(adata)
    stim = (adata.obs["condition"] == "stim").to_numpy()
    lib = np.asarray(x.sum(axis=1)).ravel()
    ngenes = np.asarray((x > 0).sum(axis=1)).ravel()
    s = species

    # the same draw as prepare.py, so eligibility here is the plan's
    drawn = prepare.draw(adata, truth)
    for cls in ("change", "no_change"):
        expect = sorted(c["gene"] for c in claims if c["comparison"] == f"{s}_lps4" and
                        c["truth"] == cls)  # fmt: skip
        assert drawn[cls] == expect, f"{s} {cls}: draw differs from claims.json"
    pool = eligible(adata, truth)
    nc = pool[truth.loc[pool, "truth"] == "no_change"]
    ch = pool[truth.loc[pool, "truth"] == "change"]
    assert (len(ch), len(nc)) == (drawn["n_change_eligible"], drawn["n_no_change_eligible"])

    # bulk: unpaired (the truth, which must reproduce bulk_truth/*.csv) and paired by index
    counts, meta = load_bulk(s)
    bulk, _ = deseq(counts, meta[["condition"]], "~condition")
    t = truth.dropna(subset=["log2FoldChange"])
    assert np.allclose(bulk.loc[t.index, "log2FoldChange"], t["log2FoldChange"], rtol=0, atol=1e-9)
    bulk_idx, _ = deseq(counts, meta[["condition", "pair"]], "~pair + condition")

    # pseudobulk on all cells, paired (Recluster's design) and unpaired
    pb = pseudobulk(x, names, key)
    pb_meta = pd.DataFrame(
        {
            "condition": [i.split("|")[1] for i in pb.index],
            "pair": [i.split("|")[0] for i in pb.index],
        },
        index=pb.index,
    )
    paired, sf = deseq(pb, pb_meta, "~pair + condition")
    unpaired, _ = deseq(pb, pb_meta[["condition"]], "~condition")
    # whole-library sizes, the denominator of CP10k, so offset = size factor minus CP10k exactly
    libsize = pd.Series(lib).groupby(key).sum()
    offset = [a - b for a, b in zip(per_animal(libsize, s), per_animal(sf, s), strict=True)]
    tables[f"{s}_offsets"] = pd.DataFrame(
        {
            "animal": [f"{s}{a}" for a in ANIMALS],
            "log2_libsize_ratio": per_animal(libsize, s),
            "log2_size_factor_ratio": per_animal(sf, s),
            "offset": offset,
        }
    )
    num[f"{s}_scaling_offset_min"] = min(offset)
    num[f"{s}_scaling_offset_max"] = max(offset)
    if s == "rat":
        for a, o in zip(ANIMALS, offset, strict=True):
            num[f"rat_scaling_offset_rat{a}"] = o

    # cycling cells
    g2m, sph = unique_present(G2M, names), unique_present(S_PHASE, names)
    num[f"{s}_n_g2m_genes"], num[f"{s}_n_s_genes"] = len(g2m), len(sph)
    cycling = (n_detected(x, names, g2m) >= MIN_MARKERS) | (
        n_detected(x, names, sph) >= MIN_MARKERS
    )
    frac = pd.Series(cycling).groupby(key).agg(["sum", "size", "mean"])
    frac.columns = ["n_cycling", "n_cells", "fraction"]
    tables[f"{s}_cycling_fraction"] = frac.rename_axis("sample").reset_index()
    if s == "rat":
        for sample, f in frac["fraction"].items():
            num[f"rat_cycling_pct_{sample.replace('|', '_')}"] = 100 * f
    num[f"{s}_cycling_pct_min"] = 100 * frac["fraction"].min()
    num[f"{s}_cycling_pct_max"] = 100 * frac["fraction"].max()

    # enrichment in cycling cells, unstimulated only, so the LPS response does not enter
    unst = ~stim
    cyc_cp = np.asarray(x[unst & cycling].sum(axis=0)).ravel() / lib[unst & cycling].sum()
    non_cp = np.asarray(x[unst & ~cycling].sum(axis=0)).ravel() / lib[unst & ~cycling].sum()
    ratio = (cyc_cp * 1e4 + PSEUDOCOUNT) / (non_cp * 1e4 + PSEUDOCOUNT)
    assoc = pd.Series(np.log2(ratio), index=names)
    assoc = assoc[~assoc.index.duplicated(keep=False)]

    # every eligible change and no-change gene: one-gene pseudobulk calls at p < 0.05
    genes = ch.union(nc).intersection(paired.index)
    assert len(genes) == len(ch) + len(nc), f"{s}: eligible genes missing from the pseudobulk"
    eg = pd.DataFrame(
        {
            "class": truth.loc[genes, "truth"],
            "bulk_log2fc": bulk.loc[genes, "log2FoldChange"],
            "bulk_lfcse": bulk.loc[genes, "lfcSE"],
            "pb_paired_log2fc": paired.loc[genes, "log2FoldChange"],
            "pb_paired_lfcse": paired.loc[genes, "lfcSE"],
            "pb_paired_p": paired.loc[genes, "pvalue"],
            "pb_unpaired_log2fc": unpaired.loc[genes, "log2FoldChange"],
            "pb_unpaired_lfcse": unpaired.loc[genes, "lfcSE"],
            "pb_unpaired_p": unpaired.loc[genes, "pvalue"],
            "cycling_assoc": assoc.reindex(genes),
        }
    ).sort_index()
    tables[f"{s}_eligible_genes"] = eg.rename_axis("gene").reset_index()
    egn, egc = eg[eg["class"] == "no_change"], eg[eg["class"] == "change"]
    num[f"{s}_n_nochange_eligible"] = len(egn)
    num[f"{s}_n_change_eligible"] = len(egc)
    for design in ("paired", "unpaired"):
        called = egn[f"pb_{design}_p"] < 0.05
        num[f"{s}_nochange_called_{design}"] = int(called.sum())
        num[f"{s}_nochange_called_{design}_pct"] = 100 * called.mean()
        same = np.sign(egc[f"pb_{design}_log2fc"]) == np.sign(egc["bulk_log2fc"])
        num[f"{s}_change_called_bulk_sign_{design}"] = int(
            ((egc[f"pb_{design}_p"] < 0.05) & same).sum()
        )
    num[f"{s}_expected_false_calls_in_10"] = 10 * num[f"{s}_nochange_called_paired"] / len(egn)
    num[f"{s}_nochange_median_lfcse_bulk"] = egn["bulk_lfcse"].median()
    num[f"{s}_nochange_median_lfcse_pb_paired"] = egn["pb_paired_lfcse"].median()
    num[f"{s}_nochange_median_lfcse_pb_unpaired"] = egn["pb_unpaired_lfcse"].median()

    ca = egn.dropna(subset=["cycling_assoc"])
    gap = ca["pb_paired_log2fc"] - ca["bulk_log2fc"]
    rho = stats.spearmanr(ca["cycling_assoc"], gap)
    num[f"{s}_spearman_cycling_vs_gap"] = float(rho.statistic)
    num[f"{s}_spearman_cycling_vs_gap_p"] = float(rho.pvalue)
    num[f"{s}_n_nochange_with_assoc"] = len(ca)
    hit = ca["pb_paired_p"] < 0.05
    top = ca["cycling_assoc"] >= 1
    num[f"{s}_nochange_2x_cycling"] = int(top.sum())
    num[f"{s}_nochange_2x_cycling_called"] = int((top & hit).sum())
    num[f"{s}_nochange_2x_cycling_called_down"] = int(
        (top & hit & (ca["pb_paired_log2fc"] < 0)).sum()
    )
    flat = ca["cycling_assoc"].abs() < 0.25
    num[f"{s}_nochange_unlinked"] = int(flat.sum())
    num[f"{s}_nochange_unlinked_called"] = int((flat & hit).sum())
    num[f"{s}_nochange_unlinked_called_up"] = int((flat & hit & (ca["pb_paired_log2fc"] > 0)).sum())
    num[f"{s}_nochange_unlinked_called_down"] = int(
        (flat & hit & (ca["pb_paired_log2fc"] < 0)).sum()
    )
    num[f"{s}_nochange_unlinked_called_median_gap"] = gap[flat & hit].median()
    num[f"{s}_nochange_unlinked_median_gap"] = gap[flat].median()

    # library sizes and genes detected per cell
    num[f"{s}_median_umi"] = float(np.median(lib))
    num[f"{s}_min_genes_per_cell"] = int(ngenes.min())
    mito = [i for i, n in enumerate(names) if n.lower().startswith("mt-")]
    mito_pct = 100 * np.asarray(x[:, mito].sum(axis=1)).ravel() / lib
    rows = []
    for cond, m in (("stim", stim), ("unst", ~stim), ("all", np.ones_like(stim))):
        rows.append(
            {
                "condition": cond,
                "n_cells": int(m.sum()),
                "median_umi": float(np.median(lib[m])),
                "median_genes": float(np.median(ngenes[m])),
                "min_genes": int(ngenes[m].min()),
                "max_mito_pct": float(mito_pct[m].max()) if mito else np.nan,
            }
        )
    tables[f"{s}_library"] = pd.DataFrame(rows)
    if s == "rat":
        num["rat_n_mito_genes"] = len(mito)
        num["rat_mito_prefixes"] = sorted({names[i][:3] for i in mito})
        num["rat_max_mito_pct"] = float(mito_pct.max())
        for r in rows[:2]:
            num[f"rat_median_umi_{r['condition']}"] = r["median_umi"]
            num[f"rat_median_genes_{r['condition']}"] = r["median_genes"]

    if s != "rat":
        return

    # the three rat false calls
    norm = pb[list(FALSE_CALLS)].div(sf, axis=0)
    cp_all = cp10k(x, lib, key, np.ones_like(stim))
    cp_non = cp10k(x, lib, key, ~cycling)
    rows = []
    for g in FALSE_CALLS:
        gi = names.get_loc(g)
        row = {
            "gene": g,
            "bulk_log2fc": bulk.at[g, "log2FoldChange"],
            "bulk_lfcse": bulk.at[g, "lfcSE"],
            "bulk_padj": bulk.at[g, "padj"],
            "bulk_by_index_log2fc": bulk_idx.at[g, "log2FoldChange"],
            "bulk_by_index_lfcse": bulk_idx.at[g, "lfcSE"],
            "bulk_by_index_p": bulk_idx.at[g, "pvalue"],
            "pb_paired_log2fc": paired.at[g, "log2FoldChange"],
            "pb_paired_p": paired.at[g, "pvalue"],
            "pb_paired_padj": paired.at[g, "padj"],
            "pb_unpaired_log2fc": unpaired.at[g, "log2FoldChange"],
            "pb_unpaired_p": unpaired.at[g, "pvalue"],
            "cycling_assoc": assoc[g],
        }
        scalings = (
            ("size_factor", norm[g]),
            ("cp10k", cp_all[gi]),
            ("cp10k_noncycling", cp_non[gi]),
        )
        for label, values in scalings:
            for a, v in zip(ANIMALS, per_animal(values, s), strict=True):
                row[f"{label}_rat{a}"] = v
        rows.append(row)
        low = g.lower()
        for k in ("bulk_lfcse", "bulk_by_index_log2fc", "bulk_by_index_lfcse", "bulk_by_index_p",
                  "pb_paired_log2fc", "pb_paired_p", "pb_paired_padj", "pb_unpaired_p",
                  "cycling_assoc"):  # fmt: skip
            num[f"rat_{low}_{k}"] = row[k]
        for label, _ in scalings:
            num[f"rat_{low}_{label}_per_animal"] = [row[f"{label}_rat{a}"] for a in ANIMALS]
    fc = pd.DataFrame(rows)
    tables["rat_false_calls"] = fc
    recluster = [provenance(s, g) for g in FALSE_CALLS]
    if all(recluster):
        ours = fc[["pb_paired_log2fc", "pb_paired_p", "pb_paired_padj"]].to_numpy()
        theirs = [[r["pb_log2fc"], r["pb_p"], r["pb_padj_genomewide"]] for r in recluster]
        assert np.allclose(ours, theirs, rtol=1e-4, atol=0), "differs from Recluster's as_labelled"
        print("rat false calls: all-cells pseudobulk equals Recluster's as_labelled")
    else:
        print("rat false calls: provenance.json not found, Recluster comparison skipped")
    num["rat_false_calls_called_unpaired"] = int(
        (unpaired.loc[list(FALSE_CALLS), "pvalue"] < 0.05).sum()
    )
    num["rat_ezh2_cycling_fold"] = float(2 ** assoc["Ezh2"])
    gi = names.get_loc("Ezh2")
    num["rat_ezh2_cycling_fold_no_pseudocount"] = float(cyc_cp[gi] / non_cp[gi])

    cc = pd.DataFrame(
        {
            "set": ["G2M" if g in g2m else "S" for g in g2m + sph],
            "bulk_log2fc": truth.loc[g2m + sph, "log2FoldChange"].to_numpy(),
            "bulk_padj": truth.loc[g2m + sph, "padj"].to_numpy(),
            "bulk_class": truth.loc[g2m + sph, "truth"].fillna("").to_numpy(),
            "pb_paired_log2fc": paired.loc[g2m + sph, "log2FoldChange"].to_numpy(),
            "pb_paired_padj": paired.loc[g2m + sph, "padj"].to_numpy(),
        },
        index=pd.Index(g2m + sph, name="gene"),
    )
    tables["rat_cell_cycle_genes"] = cc.reset_index()
    for g in CYCLE_GENES:
        low = g.lower()
        num[f"rat_{low}_bulk_log2fc"] = cc.at[g, "bulk_log2fc"]
        num[f"rat_{low}_bulk_padj"] = cc.at[g, "bulk_padj"]
        num[f"rat_{low}_bulk_class"] = cc.at[g, "bulk_class"]
        num[f"rat_{low}_pb_paired_log2fc"] = cc.at[g, "pb_paired_log2fc"]


def plain(v: object) -> object:
    if isinstance(v, list):
        return [plain(i) for i in v]
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating | float):
        return float(f"{float(v):.6g}")
    return v


def main() -> None:
    OUT.mkdir(exist_ok=True)
    claims = json.loads((HERE / "claims" / "claims.json").read_text())
    num: dict = {}
    tables: dict[str, pd.DataFrame] = {}
    for species in SPECIES:
        analyse(species, claims, num, tables)
        print(species, "done")
    for name, df in tables.items():
        df.to_csv(OUT / f"{name}.csv", index=False, float_format="%.6g")
    (OUT / "numbers.json").write_text(json.dumps({k: plain(v) for k, v in num.items()}, indent=1))
    print(json.dumps({k: plain(v) for k, v in num.items()}, indent=1))


if __name__ == "__main__":
    main()
