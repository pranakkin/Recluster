"""Download the Kang et al. 2018 PBMC dataset (38 MB) into data/ and verify its checksum.

Usage: python scripts/download_kang.py [data_dir]
"""

from __future__ import annotations

import sys
from pathlib import Path

from recluster.bench.kang import KANG_URL, check_kang, download_kang, load_kang


def main() -> None:
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data")
    path = download_kang(data_dir / "kang_2018.h5ad")
    info = check_kang(load_kang(path))
    print(f"{path}: {info['n_cells']:,} cells x {info['n_genes']:,} genes from {KANG_URL}")
    print(f"donors: {len(info['donors'])}; conditions: {info['conditions']}")
    print(f"cell types: {', '.join(sorted(info['cell_types']))}")


if __name__ == "__main__":
    main()
