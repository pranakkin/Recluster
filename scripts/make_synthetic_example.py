"""Write the synthetic example dataset used by examples/claims/synthetic_*.yaml.

Usage: python scripts/make_synthetic_example.py [out.h5ad]
"""

from __future__ import annotations

import sys
from pathlib import Path

from recluster.simulate import SimConfig, simulate


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/synthetic.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata = simulate(SimConfig(seed=0))
    truth = adata.uns["truth"]
    adata.write_h5ad(out)
    print(f"wrote {out}: {adata.n_obs} cells x {adata.n_vars} genes")
    print("Mono markers:", truth["markers"]["Mono"][:5])
    print("planted DE genes (up in Mono, stim):", truth["de_genes"])
    print("null genes:", truth["null_genes"])


if __name__ == "__main__":
    main()
