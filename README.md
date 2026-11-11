# Recluster

An independent stress test for single-cell RNA-seq claims.
For anyone who has to decide whether a single-cell result, produced by a person or by an AI agent, holds up.
Everyone is building AI that does the analysis. Recluster tells you whether to believe it.

> **Status: Milestone 1 (core engine) complete.** The benchmark on real data (Kang et al. 2018),
> the optional LLM layer, GitHub Pages reports and the full README are the next milestones.
> See [PLAN.md](PLAN.md) and [CHANGELOG.md](CHANGELOG.md).

## What it does

You give Recluster an AnnData file with raw counts and a claim:

- *"Cell population X exists."*
- *"Cell type X is more abundant in condition A than B."*
- *"Gene G is up in cell type X in condition A vs B."*

It tries to break the claim by re-running the analysis across a bounded grid of defensible
alternative choices (QC thresholds, doublet removal, normalisation, number of HVGs and PCs,
clustering resolution, batch correction, random seeds) and reports how often, and under which
choices, the claim survives. Composition and expression claims are tested across biological
replicates, never across cells. The result is one self-contained HTML report per claim with a
verdict (`ROBUST`, `FRAGILE`, `NOT SUPPORTED`), a survival chart, the choices that broke the
claim, a replicate-aware vs cell-level comparison, a methods paragraph, and a full provenance log.

Every number comes from deterministic code. No LLM is involved in any computation.

## Install

Requires Python 3.11 and [uv](https://docs.astral.sh/uv/).

```bash
git clone <this repo> && cd recluster
uv sync
```

## Quickstart (synthetic example, under two minutes)

```bash
uv run python scripts/make_synthetic_example.py
```

```bash
uv run recluster run --data data/synthetic.h5ad --claim examples/claims/synthetic_de_true.yaml
```

Open `results/<claim-name>/report.html`. The folder also holds `provenance.json`, `runs.csv`
and a copy of the claim. The other example claims in `examples/claims/` cover each claim type
and include one that is false by construction.

## Writing a claim

```yaml
claim_type: differential_expression      # population_exists | composition_shift | differential_expression
name: "ISG15 up in CD14 monocytes after IFN-beta"
columns:
  condition: stim                        # .obs column with the condition
  sample: [ind, stim]                    # replicate id; a list is combined into one key
  batch: null                            # optional; enables the batch-correction axis
  labels: cell_type                      # optional; existing labels
target:
  markers: [CD14, LYZ]                   # and/or label_value
  label_value: "CD14+ Monocytes"
comparison:
  group_a: stim
  group_b: ctrl
  direction: up_in_a
genes: [ISG15]
```

Grid axes and verdict thresholds are explicit and can be overridden with `--grid grid.yaml`
and `--thresholds thresholds.yaml`; every report prints the values that were used.

## Optional LLM layer

Everything above works with hand-written YAML and no API key. Two conveniences use the
Anthropic API when `ANTHROPIC_API_KEY` is set (in the environment or in a gitignored `.env`):

```bash
uv run recluster parse "ISG15 is up in CD14 monocytes after IFN-beta" --data data/kang_2018.h5ad --out claim.yaml
```

The model sees only the `.obs` column names and values and which gene symbols from your sentence
exist in the data. Its proposal is validated and checked against the file, shown to you labelled as
LLM output, and written only after you confirm (`--run` then runs it). No number anywhere in a
report comes from the model.

```bash
uv run recluster run --data data/kang_2018.h5ad --claim claim.yaml --summarize
```

adds a paragraph headed "LLM summary of computed results" to the report, generated from the
run's numbers. A guard discards any summary that contains a number not present in the computed
results. `recluster summarize results/<claim>/provenance.json` does the same for a finished run.
The default model is `claude-opus-5-5` (`--model` or `RECLUSTER_LLM_MODEL` to change it).

## Troubleshooting

**`No module named recluster` on macOS after `uv sync`.** If the repository lives inside
`~/Desktop` (or another folder where macOS flags dot-directories as hidden), the `.venv` files
inherit the Finder hidden flag and recent CPython releases skip hidden `.pth` files. Run:

```bash
chflags -R nohidden .venv
```

## Built by

Pranavi Akkinepally, Johns Hopkins.

MIT licensed.
