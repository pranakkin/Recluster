# Changelog

All notable changes to Recluster. Dates are when the work landed in this repository.

## [0.7.0] — unreleased

### Added
- Agent case study (`bench/agent_audit/`, harness `recluster/bench/agent_audit.py`, tests in
  `tests/test_agent_audit.py`): Claude Opus 5.5 analysed Kang 2018 through the API with one
  Python tool, and the claims in its final answer were audited with `recluster audit` under a
  protocol written before the run, with nine amendments written after reading the answer (65
  specs: 47 primary, 7 secondary, and the 11 first-run population specs, kept as superseded after
  amendment 9 re-ran them with the sample column). Results, including where the agent's
  summary disagrees with its own tables, in `bench/agent_audit/RESULTS.md`; four of the reports
  are published under `docs/reports/agent/`.
- `recluster demo`: a true claim (ISG15 up in CD14+ monocytes after IFN-β) and a claim that is
  false by construction (ISG15 differs between two fake groups of control donors), run on the
  full default grid on a seeded subset of Kang 2018 with 200 cells from each donor x condition
  sample (`recluster/demo.py`). It exits with status 1 when a claim does not get its expected
  verdict.
- `--jobs/-j` on `recluster run` and `recluster audit`: the grid runs of one claim are spread
  over worker processes, with results identical to a single process
  (`tests/test_parallel_runs.py`).
- `scripts/verify_reproduction.sh`: clones the repository into a temporary folder, installs it
  with `uv sync --frozen`, runs `recluster demo` and prints a SHA-256 hash of the key numbers
  (`scripts/demo_hash.py`, hash recorded on macOS arm64).

### Changed
- The sequential grid (`--jobs 1`) now runs with BLAS and numba limited to one thread, as the
  worker processes already did. On Linux, numpy's OpenBLAS gave slightly different results
  for different thread counts, so `--jobs 1` and `--jobs 2` disagreed in CI; on macOS they
  already matched.
- Reports name the groups in the collection-time difference ("B minus A" when group_a is B);
  before, they always printed "A minus B", which was wrong for the planted time cases, where
  group_a is the afternoon group B.
- Removed dead code (unused helpers, parameters, fields and imports) and comments that
  restated the code. The demo's key numbers hash identically before and after, and a full
  benchmark re-run gives the same verdicts and statistics (only the case-12 label above changed).
- README, HOW_IT_WORKS and bench docs: corrected stale statements (grid runs are parallel
  with `--jobs`; both held-out doublet blocks were missed; two cluster descriptions in
  RESULTS.md), added the limitation that ROBUST does not mean true, and added references.

## [0.6.0] — 2026-11-11 — Doublet populations, reproducible design measurements

### Added
- A doublet-population step in each run with doublet removal on
  (`evaluators/doublet_cluster.py: doublet_population_test`, `GridConfig.doublet_cluster_removal`,
  `GridConfig.doublet_cluster_auroc_max = 0.80`): a cluster that a cross-validated classifier
  cannot separate from simulated doublets of other cell-type groups is removed. Method and
  threshold were pre-specified in PLAN.md and chosen over three other candidates on development
  data; reports, runs.csv and provenance record removed clusters, statistics and parents.
- Two held-out planted-doublet benchmark cases (NK + CD14 Mono, CD8 T + B), defined before the
  step was built and run once: both missed (bench/RESULTS.md, "Doublet populations").
- `scripts/design_measurements.py` recomputes the numbers behind earlier design decisions and
  writes `bench/design_measurements.md` (22 rows, all agreeing with the cited values).

### Changed
- Three cited figures corrected to their reproducible values: the minimum-over-samples AUROC
  rejects about 32 % (not 24 %) at true AUROC 0.90; CD4 T Jaccard 0.90-0.97 (not 0.91-0.96);
  the planted SQSTM1 effect at half-amplitude 0.4 has padj 0.065 (not 0.07).

## [0.5.0] — 2026-11-11 — Robustness fixes from the first benchmark, and the time-of-day check

### Added
- Time-of-day confound check (`recluster/time_confound.py`) on every comparison claim, at the
  level of the sample or, when samples nest in a sample-key column whose levels each sit in one
  condition (and one collection time) and hold several samples, that level (for example the
  donor). Mode A (collection times given): permutation balance test on the groups' mean
  (cos, sin) vectors (exact up to 20,000 label assignments, otherwise 5,000 seeded draws) with
  the smallest attainable p reported; the groups count as imbalanced when the test is
  significant or the circular means differ by at least `Thresholds.time_imbalance_hours`
  (2.0 h), because with three or fewer units per group the test cannot reach alpha. Then a
  refit of the claim's replicate-level model with the same design the verdict used plus cos/sin
  covariates (Cook's filtering off, redundant time columns dropped). A condition effect that is
  not estimable after adjustment (too few replicates, or time collinear with the condition) is
  reported as inseparable from time and caps the verdict at FRAGILE, never NOT SUPPORTED;
  "held before adjustment" is the baseline run's own outcome. Mode B (times unknown): per-unit
  expression of 13 canonical core clock genes (Takahashi 2017; Wucher et al. 2023), per-gene
  tests for information and a decisive clock-state index, E-box arm minus activator arm, with a
  permutation test; CRY1 is tested gene by gene but excluded from the index because its
  transcript peaks between the arms (Ukai-Tadenuma et al. 2011); genes that differ without an
  index difference leave the label at NO TIME CONFOUND DETECTED and are listed; a test that
  cannot reach alpha makes Mode B NOT ASSESSABLE. A red-flag screen, not a phase estimate.
  Labels NO TIME CONFOUND DETECTED / POSSIBLE TIME CONFOUND / NOT ASSESSABLE, folded into the
  claim verdict as documented: the effect removed by adjustment → NOT SUPPORTED; every other
  POSSIBLE TIME CONFOUND → capped at FRAGILE.
- `Thresholds.time_imbalance_hours` (2.0 h).
- `recluster/resources/daynight_genes.csv` from Wucher et al. 2023 S8 Table (CC BY 4.0) via
  `scripts/fetch_daynight_genes.py`; claimed genes with documented day-night variation are
  flagged in the report (gene-level, not tissue-specific).
- Spec fields `columns.collection_time`, `tissue` and `target.marker_logic` (`any` | `all`).
- Paired replicate tests when every donor appears in both groups: paired t-test on logit
  proportions, `~ pair + condition` DESeq2 design (`evaluators/pairing.py`); the report labels
  the tests as paired or unpaired according to the design the run used.
- Benchmark: two planted time-of-day cases on Kang (times unknown, times known), the signal
  planted into all cells of each donor with one sample per donor, and a time column in the
  results table; the true IFN-β claims report their time-check outcome.
- Report section for the time check (balance test with its floor, refit, clock panel); the
  population panel shows the cluster-level mixture diagnostic; the methods paragraph covers the
  new rules and names the design the baseline run used.
- Determinism tests that run the engine twice and compare the full provenance apart from
  timestamps and runtimes, for a population claim and for DE claims with and without collection
  times.

### Changed
- QC thresholds are capped per dataset: a value that alone would remove more than 25 % of cells
  is pulled back to the 25 % quantile and the run is labelled "capped". Grid points that become
  identical after capping are dropped before running and recorded in provenance (`merged_runs`).
- Marker-based matching judges a cluster by a stratified within-sample AUROC (the per-sample
  Mann-Whitney U statistics summed over the summed n1 × n0; every sample with at least one
  cluster cell and one other cell contributes, with no per-sample size requirement) and requires
  more than half of the markers to be individually enriched; a cluster no sample can judge
  falls back to the whole-dataset AUROC, and the run records which statistic judged each
  candidate (`judged_by`, `marker_basis`). Clusters qualify by markers or by labels when a claim
  gives both, but by labels only when the reference comes from the claim's label column
  (`reference_is_labelled`); Jaccard is computed over cells that passed the run's QC with the
  removed fraction reported.
- Composition and expression tests require the population in at least two samples of both
  groups; otherwise the run is reported as not testable, counts as a failure, and the verdict
  states how many runs were not testable.
- Scrublet falls back to calling the top expected fraction (0.8 % per 1,000 cells of the scored
  subset) when its automatic threshold calls almost nothing.
- The cluster-level doublet mixture test is reported as a diagnostic and never gates a verdict;
  `doublet_residual_ratio_max` and `doublet_libsize_ratio_min` are informational.
- The planted-doublet benchmark rule is strict: caught only if the verdict is not ROBUST, the
  claim held in at least one run, and every run in which it held had doublet removal off.

## [0.4.0] — 2026-11-11 — Documentation, CI and GitHub Pages

### Added
- `docs/HOW_IT_WORKS.md`: plain-language explanation of every design choice and statistical test.
- `docs/index.html` plus the 16 benchmark reports under `docs/reports/`, generated by
  `scripts/build_docs.py` from `bench/results.json`; served by GitHub Pages from `/docs` with no
  build step (`.nojekyll`).
- GitHub Actions CI (`.github/workflows/ci.yml`): ruff check and format, mypy, pytest on push.
- `examples/claims/kang_isg15_cd14_monocytes.yaml` and `scripts/download_kang.py` so a single
  Kang claim can be run without the whole benchmark.
- README rewritten: what/who/why, report screenshot, Kang quickstart, benchmark headline table,
  design principles, limitations, roadmap.

## [0.3.0] — 2026-11-11 — Milestone 3: optional LLM layer

### Added
- `recluster parse "<claim>" --data X.h5ad`: an LLM (Anthropic API, key from
  `ANTHROPIC_API_KEY` or a gitignored `.env`) proposes a claim spec from plain English. The
  model sees only the `.obs` columns with their values and which gene symbols from the claim
  exist in `var_names`; its proposal is validated and checked against the AnnData, errors are
  fed back once, and the YAML is shown for confirmation (labelled as LLM output) before it is
  written or run (`--run`).
- `--summarize` on `recluster run` / `audit`, and `recluster summarize provenance.json`: a
  plain-English paragraph labelled "LLM summary of computed results", generated from a numeric
  digest of the finished run. A guard rejects any summary containing a number that is not in
  the digest, so the text can never introduce a figure the code did not compute.
- Default model `claude-opus-5-5`, overridable with `--model` or `RECLUSTER_LLM_MODEL`.
- Everything works without a key: the commands explain how to write the YAML by hand, and a
  failed summary never touches the computed outputs.

## [0.2.0] — 2026-11-11 — Milestone 2: benchmark

### Added
- `recluster bench` and the `recluster.bench` package: Kang et al. 2018 real-data claims (five
  true, one with no known effect, four false by construction via donor-level shuffles and a
  within-condition donor split), four planted errors on PBMC3k with two controls, each with a
  pre-specified pass/caught rule; results in `bench/results.json`, `bench/results_table.md`
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
