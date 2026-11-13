"""Run `recluster run` on the LPS claims drawn by prepare.py (default grid and thresholds).

Claims whose output folder already has a provenance.json are skipped, so an interrupted run can
be resumed. Writes out/summary.json.

Usage: uv run python bench/ground_truth/run_claims.py [--jobs N]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from recluster.engine import slugify
from recluster.spec import ClaimSpec

HERE = Path(__file__).parent
OUT = HERE / "out"
RUN = {"mouse_lps4", "rat_lps4"}  # PLAN.md, change 3


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, default=14)
    args = p.parse_args()
    claims = json.loads((HERE / "claims" / "claims.json").read_text())
    claims = [c for c in claims if c["comparison"] in RUN]
    OUT.mkdir(exist_ok=True)
    summary_path = OUT / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    for i, c in enumerate(claims, 1):
        key = f"{c['comparison']}/{c['gene']}"
        spec = ClaimSpec.from_yaml(HERE / c["file"])
        folder = OUT / slugify(spec.name)
        if key in summary and (folder / "provenance.json").exists():
            continue
        t0 = time.perf_counter()
        cmd = [
            sys.executable, "-m", "recluster.cli", "run",
            "--data", c["data"], "--claim", str(HERE / c["file"]),
            "--out", str(OUT), "--jobs", str(args.jobs), "--quiet",
        ]  # fmt: skip
        proc = subprocess.run(cmd, capture_output=True, text=True)
        seconds = round(time.perf_counter() - t0, 1)
        if proc.returncode != 0:
            summary[key] = {"error": proc.stderr[-2000:], "seconds": seconds}
        else:
            v = json.loads((folder / "provenance.json").read_text())["verdict"]
            summary[key] = {
                "verdict": v["label"],
                "n_holds": v["n_holds"],
                "n_runs": v["n_runs"],
                "folder": str(folder.relative_to(HERE)),
                "seconds": seconds,
            }
        summary_path.write_text(json.dumps(summary, indent=1))
        print(f"[{i}/{len(claims)}] {key} ({c['truth']}): {summary[key]}", flush=True)


if __name__ == "__main__":
    main()
