"""Print the key numbers of a `recluster demo` run and a SHA-256 hash of them.

The numbers are, for each of the two demo claims: the verdict, how many runs held, and for every
grid run its cells after QC, clusters, whether the claim held (replicate and cell level) and the
pseudobulk log2 fold change rounded to 2 decimals. Two machines that reproduce the demo print the
same hash.

Usage: python scripts/demo_hash.py results/demo
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

# recorded on macOS 26 (arm64), Python 3.11, with the uv.lock of the same commit
EXPECTED = "9a07eac39f4b4f500843f20550ea1a2c1c6fcfb47c9666188bc77c733a8f597a"


def key_numbers(demo_dir: Path) -> dict:
    out = {}
    for prov_path in sorted(demo_dir.glob("*/provenance.json")):
        prov = json.loads(prov_path.read_text())
        v = prov["verdict"]
        runs = []
        for r in sorted(prov["runs"], key=lambda r: r["index"]):
            lfc = (r.get("key_stats") or {}).get("pb_log2fc")
            runs.append(
                [
                    r["label"],
                    r["n_cells"],
                    r["n_clusters"],
                    bool(r["holds"]),
                    bool(r["holds_cell_level"]),
                    None if lfc is None else round(lfc, 2),
                ]
            )
        out[prov["claim"]["name"]] = {
            "verdict": v["label"],
            "n_holds": v["n_holds"],
            "n_runs": v["n_runs"],
            "runs": runs,
        }
    return out


def main() -> None:
    demo_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "results/demo")
    numbers = key_numbers(demo_dir)
    if len(numbers) != 2:
        raise SystemExit(f"expected 2 demo claims in {demo_dir}, found {len(numbers)}")
    for name, n in numbers.items():
        print(f"{n['verdict']:>13}  {n['n_holds']}/{n['n_runs']}  {name}")
    canonical = json.dumps(numbers, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(canonical.encode()).hexdigest()
    print(f"sha256 of the key numbers: {digest}")
    if digest == EXPECTED:
        print("matches the recorded hash")
    else:
        print(f"differs from the recorded hash {EXPECTED}")
        sys.exit(1)


if __name__ == "__main__":
    main()
