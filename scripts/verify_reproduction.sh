#!/usr/bin/env bash
# Reproduce the demo from a fresh clone and print a hash of its key numbers.
#
# Usage: scripts/verify_reproduction.sh [repository URL or path]
# Needs git and uv. Downloads Kang et al. 2018 (38 MB) and takes a few minutes.
set -euo pipefail

repo="${1:-https://github.com/pranakkin/Recluster.git}"
work="$(mktemp -d)"
echo "working in $work"

git clone --quiet "$repo" "$work/Recluster"
cd "$work/Recluster"
uv sync --frozen --quiet
uv run recluster demo --no-open --data-dir "$work/data" --out "$work/demo"
uv run python scripts/demo_hash.py "$work/demo"
