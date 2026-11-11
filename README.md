# Recluster

An independent stress test for single-cell RNA-seq claims.
For anyone who has to decide whether a single-cell result, produced by a person or by an AI agent, holds up.
Everyone is building AI that does the analysis. Recluster tells you whether to believe it.

[![ci](https://github.com/pranakkin/Recluster/actions/workflows/ci.yml/badge.svg)](https://github.com/pranakkin/Recluster/actions/workflows/ci.yml)

![A Recluster report: a FRAGILE verdict for ISG15 in monocytes on Kang 2018, with the survival strip and the survival-by-choice chart](docs/report-screenshot.png)

## What it does

You give Recluster an AnnData file with raw counts and a claim:

- *"Cell population X exists."*
- *"Cell type X is more abundant in condition A than B."*
- *"Gene G is up in cell type X in condition A vs B."*

It re-runs the analysis across a bounded grid of defensible alternative choices (QC thresholds,
doublet removal, normalisation, number of variable genes and principal components, clustering
resolution, batch correction, random seeds; 21 to 24 runs by default) and reports how often, and
under which choices, the claim survives. Composition and expression claims are tested across
biological replicates, never across cells, and the naive cell-level result is shown next to the
replicate-aware one so pseudoreplication is visible rather than hidden.

The output is one self-contained HTML report per claim: a verdict (`ROBUST`, `FRAGILE`,
`NOT SUPPORTED`), a survival chart, the choices that broke the claim, the replicate-aware versus
cell-level comparison, a deterministic methods paragraph, and a full provenance log. Every number
comes from deterministic code; a language model is only ever used, optionally, to draft the claim
file from a sentence or to write a labelled summary of numbers the code already computed.

Browse the benchmark reports: **https://pranakkin.github.io/Recluster/** (GitHub Pages from `docs/`).

## Quickstart: reproduce the Kang 2018 benchmark

Requires Python 3.11 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/pranakkin/Recluster.git && cd Recluster
```
```bash
uv sync
```
```bash
uv run recluster bench --part kang
```

That downloads the Kang et al. 2018 PBMC dataset (38 MB, checksum verified) and the two original
10x lanes with demuxlet calls from GEO (60 MB), builds ten claims (true, null and false by
construction), runs each on the full grid in parallel, and writes `bench/results.json`,
`bench/results_table.md` and a report per claim under `bench/out/`. About 25 minutes on a laptop;
`--part all` adds the planted errors on PBMC3k (another 10 minutes); `--max-runs 2` gives a
10-minute smoke version.

To run one claim instead:

```bash
uv run python scripts/download_kang.py
```
```bash
uv run recluster run --data data/kang_2018.h5ad --claim examples/claims/kang_isg15_cd14_monocytes.yaml
```

Open `results/<claim>/report.html`. A synthetic example that needs no download is in
`examples/claims/synthetic_*.yaml` (`uv run python scripts/make_synthetic_example.py` first).

## Benchmark, exactly as measured

Default grid, default thresholds, nothing tuned after seeing results. Full narrative with each
miss explained: [bench/RESULTS.md](bench/RESULTS.md).

| # | Case | Expected | Verdict | Held | Outcome |
|---|---|---|---|---|---|
| 1 | Kang: ISG15 up in CD14+ monocytes, stim vs ctrl | ROBUST | FRAGILE | 16/21 | missed |
| 2 | Kang: ISG15, IFI6, IFIT3, MX1 all up in CD14+ monocytes | ROBUST | FRAGILE | 16/21 | missed |
| 3 | Kang: ISG15 up in CD4 T cells | ROBUST | FRAGILE | 15/21 | missed |
| 4 | Kang: CD14+ monocytes exist across donors | ROBUST | FRAGILE | 16/21 | missed |
| 5 | Kang: NK cells exist across donors | ROBUST | FRAGILE | 15/21 | missed |
| 6 | Kang: monocyte proportion differs stim vs ctrl (no known shift) | NOT SUPPORTED | FRAGILE | 8/21 | false alarm |
| 7 | False: ISG15 differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 1/21 | correct |
| 8 | False: LYZ differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | correct |
| 9 | False: monocyte proportion differs between fake donor groups | NOT SUPPORTED | NOT SUPPORTED | 0/21 | correct |
| 10 | False: ISG15 up in fake group A (donor-level shuffle) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | correct |
| 11 | Planted: 150 synthetic doublets as a fake population | not ROBUST | ROBUST | 23/23 | missed |
| 12 | Planted: over-clustering shard of CD4 T cells | not ROBUST | NOT SUPPORTED | 0/23 | caught |
| 13 | Control: a shard that is a whole population (B cells) | ROBUST | ROBUST | 23/23 | control OK |
| 14 | Planted: QC threshold deletes NK cells | failure table names max_mito_pct | NOT SUPPORTED | 4/23 | caught |
| 15 | Planted: batch artifact claimed as a CD4 T subtype | not ROBUST | NOT SUPPORTED | 4/24 | caught |
| 16 | Control: the real CD4 T population under the batch effect | ROBUST | ROBUST | 24/24 | control OK |

Three things to read off this table. The replicate-aware tests do their job: every false claim is
NOT SUPPORTED while the naive cell-level tests would have called three of the four in 81 to 100 %
of runs. The grid's `min_genes=500` setting is not a defensible alternative on a dataset whose
median cell has 519 genes; it alone makes every true Kang claim FRAGILE, and the reports attribute
it correctly. And the planted doublet case is a clean miss, for reasons spelled out in RESULTS.md.
Scrublet against demuxlet on the real lanes: AUROC 0.87, recall 0.27, nothing called in one lane.

## Writing a claim

```yaml
claim_type: differential_expression      # population_exists | composition_shift | differential_expression
name: "ISG15 up in CD14+ monocytes, stim vs ctrl"
columns:
  condition: label                       # .obs column with the condition
  sample: [replicate, label]             # biological replicate; a list is combined into one key
  batch: null                            # optional; enables the batch-correction axis
  labels: cell_type                      # optional; existing labels become the reference set
target:
  markers: [CD14, LYZ, S100A8]           # and/or label_value
  label_value: "CD14+ Monocytes"
comparison:
  group_a: stim
  group_b: ctrl
  direction: up_in_a
genes: [ISG15]
```

`recluster run` re-derives everything from raw counts. `recluster audit` takes an object another
tool produced (processed matrix plus raw counts in `layers['counts']` or `.raw`) and checks its
labels against a fresh derivation. Grid axes and thresholds can be overridden with `--grid` and
`--thresholds` YAML files; every report prints the values used.

## Optional LLM layer

Everything above works with hand-written YAML and no API key. With `ANTHROPIC_API_KEY` set (in
the environment or a gitignored `.env`):

```bash
uv run recluster parse "ISG15 is up in CD14 monocytes after IFN-beta" --data data/kang_2018.h5ad --out claim.yaml
```

proposes a claim file. The model sees only the `.obs` column names and values and which gene
symbols from the sentence exist in the data; the proposal is validated against the file and shown
to you, labelled as LLM output, before it is written. `--summarize` on a run adds a paragraph
headed "LLM summary of computed results"; a guard discards any summary containing a number the
code did not compute.

## Design principles

- **Deterministic stress test.** Every stochastic step is seeded and neighbour search is exact;
  re-running a config reproduces the numbers, and a test checks it.
- **Replicate-aware statistics.** Donors, not cells, are the unit of inference (Welch t-test on
  logit proportions; DESeq2 on per-donor pseudobulk). The cell-level result is shown for contrast
  and raises a pseudoreplication flag when it disagrees.
- **Provenance.** Package versions, input hash, grid, thresholds and every run's parameters and
  outputs are written to `provenance.json` and embedded in the report.
- **LLM only at the edges.** It may draft a claim file (checked against the data, confirmed by
  you) and summarise computed numbers (guarded). It never touches a computation.
- **Explicit thresholds.** All cutoffs are configurable and printed; changes are logged with
  reasons in `PLAN.md`.

[docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md) explains every design choice and test.

## Limitations

- The QC axes are absolute counts, not adapted to the dataset; on shallow data one setting can
  dominate the verdict (it did, on Kang).
- Existence is scored with a Jaccard that penalises cells removed by QC, so "QC deleted part of
  it" and "it is not there" are only distinguished in the failure table, not in the verdict.
- Populations are matched by marker enrichment across both conditions at once; a treatment
  that remodels marker expression can make a population drop out of one condition and
  fabricate a composition shift.
- Multi-marker claims use an additive marker score, which cannot express co-expression.
- Doublet removal relies on Scrublet's automatic threshold, which failed on one real lane and on
  a block of identical synthetic doublets.
- Composition tests are unpaired even for paired designs; runs within a claim are sequential.
- The LLM layer's two live API calls were built against the current SDK reference and tested
  with fake clients, not against the live API.

## Roadmap

The main direction is **audit mode for agent outputs**: take the object an analysis agent
hands back, re-derive from its stored counts, and report which of its labels, abundances and
marker calls survive. Concretely, next:

1. Label-seeded and per-condition population matching, and a "present in both groups"
   precondition for composition and expression tests.
2. Data-relative QC axes (quantiles of genes per cell) and existence scored on the cells a run
   kept, with the QC-removed fraction reported separately.
3. A co-expression marker score and a heterotypic-marker doublet check; a Scrublet threshold
   fallback tied to the expected doublet rate.
4. Paired replicate tests when the design is paired; grid runs in parallel.
5. An audit summary across all labels in an object, not one claim at a time.

## Repository layout

```
recluster/            engine: spec, pipeline, grid, matching, evaluators, verdict, report, llm, bench
bench/                benchmark runner outputs: RESULTS.md, results.json, results_table.md
docs/                 GitHub Pages: index.html, reports/, HOW_IT_WORKS.md, report-screenshot.png
examples/claims/      claim YAML examples (synthetic, PBMC3k, Kang)
scripts/              make_synthetic_example.py, make_processed_example.py, download_kang.py, build_docs.py
tests/                70 offline tests on synthetic data
```

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
