# Changelog

All notable changes to Recluster. Dates are when the work landed in this repository.

## [0.2.0] — 2026-11-11 — Milestone 2: benchmark

### Added
- `recluster bench` and the `recluster.bench` package: Kang et al. 2018 real-data claims (five
  true, one with no known effect, four false by construction via donor-level shuffles and a
  within-condition donor split), four planted errors on PBMC3k with two controls, each with a
  pre-registered pass/caught rule; results in `bench/results.json`, `bench/results_table.md`
  and the narrative `bench/RESULTS.md`. Claims run in parallel worker processes (`--jobs`).
- Doublet detection benchmark against demuxlet calls on the two original Kang 10x lanes from
  GEO (GSE96583 batch 2), with per-lane precision, recall and AUROC.
- Reports show the matched population's median Scrublet score next to the overall median
  (informational; not part of the verdict).

### Changed
- Scrublet is run per sample by Recluster's own loop with exact neighbours; scanpy's
  `batch_key` path fails with scipy 1.17 on large samples and was also slower.
- The `max_mito_pct` axis is dropped from the grid when the input has no mitochondrial genes
  (it would only duplicate runs); the methods paragraph says so.
- `--max-runs` is accepted by `recluster bench` for quick reproductions.

## [0.1.0] — 2026-11-10 — Milestone 1: core engine

### Added
- Claim spec (`recluster/spec.py`): pydantic models for `population_exists`,
  `composition_shift` and `differential_expression`, YAML loading, and validation against
  the AnnData (columns, label values, condition groups, genes, replicate nesting).
- Deterministic pipeline (`recluster/pipeline.py`): QC → Scrublet (per sample) → log-normalise
  or Pearson residuals → HVG → PCA → optional Harmony → exact kNN → Leiden (igraph flavour).
  QC + doublet stage cached per QC setting.
- Perturbation grid (`recluster/grid.py`): baseline + one-at-a-time + seeded random
  combinations over 9 axes; 24 runs by default; configurable from YAML.
- Population matching (`recluster/matching.py`): marker AUROC per cluster, label-precision
  matching, Jaccard against a reference set counting QC-removed cells as absent, sample spread.
- Evaluators: population recovery with a multi-sample check; composition shift with a Welch
  t-test on logit proportions and Mann-Whitney U across samples (plus the naive cell-level
  chi-square for contrast); differential expression with cell-level Wilcoxon and donor-level
  pseudobulk DESeq2 (pydeseq2), with per-gene agreement and a pseudoreplication flag.
- Verdicts with explicit, configurable thresholds (`Thresholds`), per-axis failure rates and
  one-at-a-time flips ("what broke it").
- Provenance JSON per claim: package versions, platform, input SHA-256, grid, thresholds,
  every run's parameters and outputs, runtime.
- Self-contained HTML report with inline SVG survival strip and survival-by-choice chart,
  replicate-aware vs cell-level panel, "as labelled" panel, deterministic methods paragraph.
- CLI: `recluster run`, `recluster audit`, `recluster version`.
- Synthetic data generator (`recluster/simulate.py`) with planted truth, used by the tests and
  the example claims in `examples/claims/`.
- 44 pytest tests: spec validation, each evaluator (true and false claims), grid and verdict
  rules, pseudobulk aggregation, end-to-end engine run, determinism.
