"""Write a "processed by another tool" version of the synthetic example, for `recluster audit`.

The object has log-normalised values in .X, raw counts in layers['counts'], and cluster
labels in .obs['tool_label'] produced by a plain scanpy pipeline, as an agent or a notebook
might hand over.

Usage: python scripts/make_processed_example.py [in.h5ad] [out.h5ad]
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import anndata as ad
import scanpy as sc


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/synthetic.h5ad")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/synthetic_processed.h5ad")
    a = ad.read_h5ad(src)
    a.layers["counts"] = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=2000)
        sc.pp.pca(a, n_comps=30, mask_var="highly_variable", svd_solver="arpack")
        sc.pp.neighbors(a, n_neighbors=15, transformer="sklearn")
        sc.tl.leiden(a, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    # name clusters by their majority true type, the way an annotator would
    majority = a.obs.groupby("leiden", observed=True)["cell_type"].agg(lambda s: s.mode().iloc[0])
    a.obs["tool_label"] = a.obs["leiden"].map(majority).astype(str)
    del a.obs["cell_type"]  # the audit must not see the truth
    a.write_h5ad(out)
    print(f"wrote {out}: {a.n_obs} cells, labels: {sorted(a.obs['tool_label'].unique())}")


if __name__ == "__main__":
    main()
