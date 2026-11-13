"""Run `recluster audit` on every spec from make_specs.py against the agent's annotated.h5ad.

Each claim goes through the command-line tool exactly as a user would run it. Claims whose
output folder already has a provenance.json are skipped, so an interrupted run can be resumed.
Writes out/summary.json.

Usage: uv run python bench/agent_audit/run_audits.py [--jobs N]
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
DATA = HERE / "run" / "workspace" / "annotated.h5ad"
OUT = HERE / "out"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, default=14)
    args = p.parse_args()
    items = json.loads((HERE / "claims" / "claims.json").read_text())
    OUT.mkdir(exist_ok=True)
    summary_path = OUT / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    for i, c in enumerate(items, 1):
        spec = ClaimSpec.from_yaml(HERE / c["file"])
        folder = OUT / slugify(spec.name)
        if (folder / "provenance.json").exists() and c["id"] in summary:
            continue
        t0 = time.perf_counter()
        cmd = [
            sys.executable, "-m", "recluster.cli", "audit",
            "--data", str(DATA), "--claim", str(HERE / c["file"]),
            "--out", str(OUT), "--jobs", str(args.jobs), "--quiet",
        ]  # fmt: skip
        proc = subprocess.run(cmd, capture_output=True, text=True)
        seconds = round(time.perf_counter() - t0, 1)
        if proc.returncode != 0:
            summary[c["id"]] = {"error": proc.stderr[-2000:], "seconds": seconds}
        else:
            v = json.loads((folder / "provenance.json").read_text())["verdict"]
            summary[c["id"]] = {
                "verdict": v["label"],
                "n_holds": v["n_holds"],
                "n_runs": v["n_runs"],
                "n_not_testable": v.get("n_not_testable"),
                "folder": str(folder.relative_to(HERE)),
                "seconds": seconds,
            }
        summary_path.write_text(json.dumps(summary, indent=1))
        print(f"[{i}/{len(items)}] {c['id']}: {summary[c['id']]}", flush=True)


if __name__ == "__main__":
    main()
