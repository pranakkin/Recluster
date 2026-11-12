# Recluster design log

Decisions and deviations are recorded here as they were made. The doublet-population step of
0.6.0, its adoption rules and its two held-out cases were written here and committed before the
step was built: 6201b82, 8084ad2 and 1c1e48e (the threshold) come before f1dd9ea in the git
history.

## What we are building

An independent stress test for single-cell RNA-seq claims. Input: an AnnData file plus a
structured claim (YAML). Output: one self-contained HTML report per claim with a verdict
(`ROBUST` / `FRAGILE` / `NOT SUPPORTED`), a survival chart across a grid of defensible
analysis choices, the choices that broke the claim, a replicate-aware vs cell-level
comparison, a deterministic methods paragraph, and a full provenance log.

Hard rule: every number comes from deterministic code. An LLM may only (a) draft a claim
spec from plain English and (b) summarize numbers the code already computed. Both are
labelled as LLM output.

## Environment

- Python 3.11 (uv-managed interpreter `cpython-3.11.16`), `uv` for the venv and lockfile.
- Machine used for development: Apple M5 Pro, 15 cores, 24 GB, arm64. Timings in this
  repo were measured here unless stated otherwise.
- All versions pinned in `pyproject.toml`; `uv.lock` committed.

## Module layout

```
recluster/
  __init__.py          version string
  cli.py               typer CLI: run, audit, bench (M2), parse (M3)
  spec.py              pydantic claim spec (3 claim types) + YAML loader + AnnData column check
  config.py            GridConfig (axes, baseline, design) and Thresholds (verdict cutoffs), YAML-loadable
  data.py              load AnnData, locate raw counts (.X / .raw / layers['counts']), sha256, obs validation
  pipeline.py          one deterministic pipeline run: QC -> doublets -> norm -> HVG -> PCA -> (harmony) -> kNN -> leiden
  grid.py              bounded design: baseline + one-at-a-time + seeded random combinations
  matching.py          find the target population in a run (marker AUROC, cluster matching, Jaccard)
  evaluators/
    population.py      population_exists
    composition.py     composition_shift (per-sample proportions, logit t-test + Mann-Whitney; naive chi-square for contrast)
    de.py              differential_expression (cell-level Wilcoxon vs donor-level pseudobulk DESeq2)
    pseudobulk.py      sum raw counts per sample within a cell set
  verdict.py           survival fractions -> verdict with explicit, configurable thresholds; "what broke it" table
  provenance.py        JSON log: params, package versions, seeds, input hash, runtime
  report/
    html.py            self-contained HTML (inline CSS, inline SVG chart; no JS, no external assets)
    methods.py         deterministic methods paragraph from a template
    templates/report.html.j2
  engine.py            orchestrates load -> grid -> run -> evaluate -> verdict -> report
  simulate.py          synthetic count generator (donors, conditions, planted effects) for tests and bench
  llm/                 M3: parse.py (claim -> spec YAML), summarize.py (numbers -> prose)
bench/                 M2: Kang 2018 real-data benchmark + planted errors on pbmc3k, RESULTS.md
tests/                 pytest, offline, synthetic data only
docs/                  GitHub Pages: index.html + benchmark reports + HOW_IT_WORKS.md
```

## Claim spec (YAML)

```yaml
claim_type: differential_expression        # population_exists | composition_shift | differential_expression
name: "ISG15 up in CD14 monocytes after IFN-beta"
columns:
  condition: stim                          # required for composition_shift / differential_expression
  sample: [ind, stim]                      # string or list; a list is joined into one replicate key
  batch: null                              # optional; enables the batch-correction axis
  labels: cell_type                        # optional; existing labels (required in audit mode)
target:                                    # how to find the population
  markers: [CD14, LYZ]                     # optional
  negative_markers: []                     # optional
  label_value: "CD14+ Monocytes"           # optional; value in columns.labels
  min_cells: 20
comparison:                                # composition_shift / differential_expression only
  group_a: stim
  group_b: ctrl
  direction: up_in_a                       # population: n/a; composition: higher_in_a|lower_in_a|different; DE: up_in_a|down_in_a|either
genes: [ISG15, IFI6]                       # differential_expression only
require: all                               # all | any (multi-gene DE)
```

Validation: pydantic v2 models. `spec.check_against(adata)` verifies every named column and
value exists and that each sample maps to exactly one condition (replicates must be nested
in condition; otherwise the user is told to use a list key such as `[donor, condition]`).

## Perturbation grid

Each run is a `PipelineParams`:

| axis | options (default) | baseline |
|---|---|---|
| `min_genes` | 100, 200, 500 | 200 |
| `max_mito_pct` | 5, 10, 20 | 10 |
| `doublet_removal` | off, on (scanpy's scrublet, per sample when a sample column exists) | on |
| `normalization` | `log_total` (CP10k + log1p), `pearson` (analytic Pearson residuals) | `log_total` |
| `n_hvgs` | 1000, 2000, 3000 | 2000 |
| `n_pcs` | 20, 30, 50 | 30 |
| `resolution` | 0.3, 0.5, 1.0, 1.5, 2.0 (leiden, igraph flavour) | 1.0 |
| `batch_correction` | off, on (harmonypy on PCs) — only when a batch column exists | on if batch column |
| `seed` | 0, 1, 2 | 0 |

Design: baseline + one-at-a-time (every non-baseline value of every axis) + `n_random`
seeded random combinations (default 6). Default grid is 24 runs without batch axis, 25 with.
All of this is in `GridConfig` and overridable from YAML (`--grid`).

Caching: QC + doublet scoring depends only on (`min_genes`, `max_mito_pct`), so that stage is
computed once per distinct QC setting and reused. Everything downstream is recomputed.

What the normalization axis does and does not change: it changes the representation used
for HVG selection, PCA and clustering. Marker scoring and cell-level Wilcoxon always use
CP10k log1p expression; pseudobulk always uses raw counts. This keeps the tests comparable
across runs and is stated in the methods paragraph.

## Finding the target population in a run (`matching.py`)

- Per-cell marker score: mean CP10k-log1p expression of `markers` minus mean of
  `negative_markers`, z-scored across cells.
- Per-cluster marker AUROC: AUROC of the marker score for cells in the cluster vs all others.
- Marker mode: candidate clusters = clusters with AUROC >= `marker_auroc_min` (0.80). The
  matched set M is their union. This tolerates a population splitting into sub-clusters; a
  population merged into a larger cluster lowers the AUROC and fails, which is the intent
  ("exists as a distinct population").
- Label mode (no markers): candidate clusters = clusters where >= 50 % of cells carry the
  reference label. M = union.
- Reference set R: labelled cells if `label_value` is given; otherwise M from the baseline run
  (if the baseline fails to find it, the first run that finds it; recorded in provenance).
- Jaccard(M, R) is computed over barcodes. Cells removed by QC in a run count as absent from M,
  so a QC setting that deletes the population lowers Jaccard.
- Sample spread: number of samples contributing >= `min_cells_per_sample` (5) cells to M, and
  the largest single-sample share of M.

## Evaluators and per-run `holds`

**population_exists**
- `recovered` = |M| >= `min_cells` and Jaccard >= `jaccard_min` (0.50) and (if markers) marker AUROC of M >= 0.80.
- `replicate_ok` = samples contributing >= 5 cells >= `min_samples_present` (2) when a sample column exists; else true.
- `holds` = recovered and replicate_ok.

**composition_shift**
- Per sample: proportion of that sample's cells that are in M.
- Replicate-aware test: Welch t-test on logit((k+0.5)/(n+1)) across samples, plus Mann-Whitney U on raw proportions. Effect: difference in mean log-odds and ratio of mean proportions.
- `holds` = replicate-aware p < `alpha` (0.05) and direction matches (`different` ignores direction).
- For contrast only: naive cell-level chi-square on the 2x2 table of cells. Never used for the verdict.

**differential_expression**
- Cells in M, group_a vs group_b.
- Cell-level: Wilcoxon rank-sum on CP10k-log1p per gene, BH-adjusted over the claimed genes; log2 fold change of mean normalised expression with pseudocount.
- Replicate-aware: pseudobulk by summing raw counts per sample within M; pydeseq2 with design `~ condition` on genes passing a minimum count filter; report log2FC and BH padj for the claimed genes.
- `holds` = pseudobulk padj < alpha and direction matches, combined over genes with `require`.
- `holds_cell_level` recorded for the comparison panel. Pseudoreplication flag when cell-level survival >= 0.8 and replicate-aware survival < 0.5 (Squair et al. 2021, Nat Commun 12:5692).

## Verdict (`verdict.py`)

- survival = fraction of grid runs where `holds`.
- `ROBUST` if survival >= `robust_min_survival` (0.80).
- `FRAGILE` if `fragile_min_survival` (0.20) <= survival < 0.80. The report lists the axis values whose runs failed, with per-axis-value failure rates, and the one-at-a-time runs that flipped the result relative to baseline.
- `NOT SUPPORTED` if survival < 0.20.
- Thresholds live in `Thresholds` (pydantic), are printed in every report, and can be overridden with `--thresholds thresholds.yaml`.

## Provenance

Per claim: `provenance.json` with recluster version, timestamp, Python and package versions,
platform, input path + sha256 + shape, the claim spec, grid config, thresholds, and for every
run the params, seed, cell and cluster counts, wall time and evaluator outputs. Determinism:
fixed seeds for scrublet, PCA (arpack), harmony, leiden; reported numbers rounded to 6
significant digits; a test asserts two runs produce identical results.

## Build order (Milestone 1)

1. pyproject + uv venv; confirm scanpy, igraph, harmonypy, pydeseq2 import on arm64.
2. `spec.py`, `config.py`.
3. `data.py`, `simulate.py`.
4. `pipeline.py`; smoke test on synthetic data and pbmc3k.
5. `grid.py`.
6. `matching.py`, `evaluators/*`.
7. `verdict.py`, `provenance.py`.
8. `report/*`.
9. `engine.py`, `cli.py`.
10. tests, end-to-end run, timing.

## Milestone 1 status (2026-11-10): done

- 44 tests pass (`uv run pytest`, ~45 s). ruff and mypy clean.
- Full 24-run grid on the synthetic example (1,800 cells): ~30 s for population/composition
  claims, ~45 s for DE claims (pydeseq2 adds ~1 s per run). Verdicts on the four synthetic
  example claims match their construction: population ROBUST 24/24, composition ROBUST 24/24,
  planted DE ROBUST 24/24, null DE NOT SUPPORTED 0/24.
- Determinism verified by `tests/test_engine.py::test_determinism_same_config_same_numbers`
  (two full runs produce identical provenance apart from timestamps and runtimes).
- PBMC3k (2,700 cells, no replicates, 23-run grid, ~23 s each): CD14+ monocytes ROBUST 20/23
  (fails only when resolution 0.3 merges them with FCGR3A+ monocytes); FCGR3A+ monocytes
  FRAGILE 18/23 (broken by resolution 0.3/0.5 and 1000 HVGs); megakaryocytes (15 cells)
  ROBUST 21/23, deleted by `min_genes=500`. These were not tuned; they are the first runs.
- Not yet measured: runtime on a ~25k-cell dataset (Kang). Expect roughly 10x pbmc3k for
  population/composition claims; DESeq2 on 16 pseudobulk samples x ~15k genes per run is the
  unknown. If it exceeds the 10-minute budget, run grid points in parallel (joblib) before
  reducing the grid.

## Decisions and deviations log

- 2026-11-10: Python 3.11 is not on PATH; using uv's managed 3.11.16.
- 2026-11-10: Clustering uses scanpy's leiden with `flavor="igraph"` (faster, deterministic with a seed) rather than leidenalg.
- 2026-11-10: Doublet detection uses `scanpy.pp.scrublet` (scanpy's native implementation of Scrublet) so no separate `scrublet` package is needed. Its automatic threshold needs `scikit-image`, which is therefore a pinned dependency.
- 2026-11-10: Scrublet is run once per QC setting with a fixed seed (`scrublet_seed`), so the `doublet_removal` axis is purely on/off and the `seed` axis covers PCA, Harmony and Leiden only.
- 2026-11-10: kNN uses the exact sklearn transformer rather than pynndescent, for determinism on any machine. Costs a few seconds on 25k cells.
- 2026-11-10: The baseline turns batch correction on whenever a batch column is declared; `build_grid` forces it off (and drops the axis) when there is none.
- 2026-11-10: DE multiplicity is adjusted over the claimed genes (a claim names genes a priori). The genome-wide BH value from DESeq2 is also reported. This choice is stated in the methods paragraph.
- 2026-11-10: Composition verdict uses the Welch t-test on logit proportions by default (`Thresholds.composition_test`); Mann-Whitney is reported alongside and can be made decisive.
- 2026-11-10: Reports use inline SVG and inline CSS only; `tests/test_engine.py` asserts no scripts, stylesheets or remote URLs.
- 2026-11-10: macOS quirk: a `.venv` created inside `~/Desktop` gets the Finder "hidden" flag on every file, and CPython ≥ 3.11.16 / 3.12.x skips hidden `.pth` files, which breaks the editable install (`No module named recluster`). Fix: `chflags -R nohidden .venv`. Documented in README troubleshooting. Not a uv bug (a venv created elsewhere has no flag).

## Milestone 2 status (2026-11-11)

- `recluster bench` reproduces everything in ~35 min (4 workers). Results in `bench/RESULTS.md`.
- Score against pre-specified rules: 4/4 false claims correct (cell-level stats called 3 of
  them in 81-100 % of runs), 3/4 planted errors caught, 2/2 controls OK, 5/5 true Kang claims
  FRAGILE instead of ROBUST (all because `min_genes=500` removes 44 % of this shallow dataset),
  1 planted error missed (block of identical doublets evades Scrublet), 1 false alarm (null
  composition claim becomes FRAGILE under the same QC setting).
- Candidate v0.3 fixes, not applied: data-relative QC axes; existence scored on kept cells;
  Scrublet threshold from expected doublet rate; heterotypic-marker doublet check.

### Milestone 2 decisions (2026-11-11)

- Kang 2018 source: the pertpy-distributed object (`exampledata.scverse.org/pertpy/kang_2018.h5ad`, 38 MB, SHA-256 pinned in `recluster/bench/kang.py`), downloaded directly rather than installing pertpy. It has raw counts, `label`, `replicate`, `cell_type`. It has **no mitochondrial genes** (filtered upstream) and **no demuxlet multiplet column** (singlets only).
- Because the object has no mito genes, `build_grid` now drops the `max_mito_pct` axis when the input has no mito genes (it would only duplicate runs). Kang therefore runs a 21-run grid (14 one-at-a-time + 6 random), PBMC3k 23 runs. The methods paragraph states this.
- Doublet check: done against demuxlet on the two original 10x lanes from GEO (GSE96583 batch 2, six files, SHA-256 pinned in `recluster/bench/doublets.py`), not on the processed object. demuxlet only sees inter-individual doublets, so agreement is a lower bound.
- scanpy's Scrublet failed with scipy 1.17 on large subsets (`csr_array has no attribute getnnz`, from its approximate-neighbour path). Scrublet is now run per sample by Recluster's own loop with `use_approx_neighbors=False`. Exact neighbours are deterministic and, on Kang, six times faster.
- Probing the planted-doublet construction before the benchmark showed that 150 near-identical mono+B doublets evade Scrublet (they are each other's neighbours; median score 0.10 vs 0.04). Decided before the run: reports now show the matched population's median Scrublet score against all cells as an informational number. It does not enter the verdict. Thresholds unchanged.
- PBMC3k reference labels for building planted errors use resolution 1.0 (0.5 merged CD8 T into NK).
- Benchmark claims run in parallel worker processes (`--jobs`, default 4); the grid itself is still sequential within a claim. (Since 41d27d9, `recluster run` and `recluster audit` also take `--jobs`, which spreads one claim's grid runs over processes; `bench --jobs` still sets how many claims run in parallel.)

## v0.5 plan (2026-11-11): robustness fixes and Milestone 2.5 (time-of-day confound)

Robustness, each change motivated by a measured v0.1 failure mode (benchmark re-run, v0.1 table kept). Where the design moved while building, the bullet says what it was first and why it changed:
- `GridConfig.qc_max_cell_loss`: a QC threshold that alone would remove more than 25 % of cells is pulled back to the quantile that removes exactly 25 %; the run records both values. Grid points that become identical after capping are dropped before running and recorded in provenance (`merged_runs`). Reason: two requested thresholds that cap to the same value would run the same pipeline twice and count one configuration twice in the survival rate.
- Matching: a cluster's marker AUROC is a stratified within-sample AUROC (the Mann-Whitney U of every sample holding at least one cluster cell and one other cell, summed, over the summed n1 × n0); a cluster no sample can judge, or every cluster without a sample column, is judged by the whole-dataset AUROC, and the statistic that judged each candidate is recorded (`judged_by`, `marker_basis`). Reason: the first implementation took the minimum over samples holding ≥ 10 cluster cells; on simulations it rejected a population with true AUROC 0.85 and 12 cells in each of 16 samples 95 % of the time (stratified statistic: 0 %, mean 0.851) and still rejected about 32 % at true AUROC 0.90.
- Matching: a cluster qualifies by markers or by labels (`Thresholds.candidate_rule`); by markers only when more than half of the markers are individually enriched (per-gene AUROC ≥ 0.6); `marker_logic: all` for co-expression claims (minimum score); Jaccard over cells that passed the run's QC with the QC-removed fraction reported. Reason: v0.1's whole-dataset AUROC moved with the cell-type mix a QC setting kept, an additive score could not express co-expression, and QC deletion was indistinguishable from absence.
- `reference_is_labelled` (`match_population`): a cluster can qualify by labels only when the reference comes from the claim's label column; a reference carried over from the baseline run's match is used only for the Jaccard overlap. Reason: without it a merged low-resolution cluster holding most of the baseline population qualified "by labels" with no marker evidence.
- Cluster-level doublet check (`evaluators/doublet_cluster.py`): mixture-of-two-clusters residual and library-size ratio, shown in the report's population panel as a diagnostic; `doublet_residual_ratio_max` and `doublet_libsize_ratio_min` only label it and no verdict depends on them. Reason: the plan had it gate a verdict; no cutoff tried separated true mixtures from genuine intermediate populations with a margin worth enforcing.
- Scrublet threshold fallback at the expected doublet rate (0.8 % per 1,000 cells of the scored subset) when the automatic threshold calls almost nothing. The expectation equals the 10x loading rule only when each sample had its own lane; measured on the two Kang GEO lanes with sample = lane it is 10.8-11.1 % against 11.5 % called by demuxlet. Reason: the automatic threshold called nothing on one Kang lane.
- Paired replicate tests (`evaluators/pairing.py`): paired t-test and `~ pair + condition` DESeq2 design when every donor appears in both groups. Measured on Kang (8 labelled types, 7 testable): the unpaired Welch test calls none (all p ≥ 0.67); the paired test calls CD14+ monocytes (p = 0.035, 23.8 % vs 25.4 %, lower in stim for 7 of 8 donors; the benchmark's null case, reported as measured); standard errors shrink 2-9 fold. Reason: Kang is a split-sample design and the unpaired test ignores the pairing.
- Composition and expression tests require the population in ≥ 2 samples of both groups; otherwise the run is not testable, counts as a failure, and the verdict says how many runs were. Kept as designed. Reason: when the matching finds the population in one condition only, a biological absence cannot be told from a clustering that lost that condition's version of the population, and v0.1's only composition false alarm (a 0.1 % vs 35 % "shift") arose exactly so; the cost, that a condition-specific population cannot be tested as a composition shift, is documented as a limitation.
- Planted-doublet benchmark rule: caught only if the verdict is not ROBUST, the claim held in at least one run, and every holding run had doublet removal off. Reason: a claim that never holds says nothing about doublet removal, so a "never holds counts as caught" reading is not allowed.

Milestone 2.5 (`recluster/time_confound.py`):
- Spec: `columns.collection_time` (hours 0-24 or "HH:MM"), `tissue` (optional; recorded, but the public list is gene-level and cannot use it).
- Units: the sample, or the coarsest sample-key column whose levels each sit in one condition (and one collection time) and hold several samples, with ≥ 2 levels per group; `sample = [donor, label]` with a condition assigned per donor makes the donor the unit. Reason: the plan tested "across samples"; two samples of one donor are not two independent observations of when the donor was sampled.
- Mode A balance: permutation test on the distance between the groups' mean (cos, sin) vectors (exact when ≤ 20,000 label assignments, otherwise 5,000 seeded draws), the smallest attainable p reported (2/C(n,k) for equal groups, 1/C(n,k) otherwise); imbalanced when p < alpha or the circular means differ by ≥ `Thresholds.time_imbalance_hours` (2.0 h). Reason: with three or fewer units per group the test cannot reach 0.05 (3 vs 3: floor 0.10), and a 2 h shift moves a rhythmic gene by up to 2·sin(π·2/24) of its half-amplitude, 26 % of its peak-to-trough swing, which should not pass as balanced.
- Mode A adjustment: the same design the verdict used (`~ condition` or `~ pair + condition`) plus cos/sin, Cook's filtering off so every claimed gene keeps a p in both fits, redundant time columns dropped (two distinct times make sin an affine function of cos; one time per pair is absorbed by the pairing); a condition coefficient that is not estimable (too few replicates for the terms, or time collinear with the condition, e.g. every control sample at one hour and every treated sample at another) is reported as "cannot be separated from time" and caps the verdict at FRAGILE, never NOT SUPPORTED. Reason: the plan reported shrinkage only; a collinear design has no shrinkage to report and must not read as a vanished effect.
- Mode A decision: balanced → NO TIME CONFOUND DETECTED; imbalanced and the claim did not hold anyway → POSSIBLE (cap); not estimable → POSSIBLE (cap); held, adjusted p ≥ alpha and (shrinkage ≥ `time_shrink_flag` 50 % or a flip) → POSSIBLE with the effect removed → NOT SUPPORTED; adjusted p ≥ alpha with the estimate barely moved → POSSIBLE (cap; degrees of freedom spent on time, cannot separate); adjusted significant but shrunk ≥ 50 % → POSSIBLE (cap); otherwise NO TIME CONFOUND DETECTED. "Held before adjustment" is the baseline run's own outcome. Reason: the refit must never contradict the grid.
- Mode B: per-unit log2 CPM of the 13 core clock genes (Takahashi 2017; the panel Wucher et al. 2023 report as day-night variable); per-gene Welch t with BH for information; decisive clock-state index = mean z of the E-box arm (PER1, PER2, PER3, CRY2, NR1D1, NR1D2, DBP, TEF, CIART) minus mean z of the activator arm (ARNTL, CLOCK, NPAS2), permutation test on the mean difference; CRY1 stays in the per-gene panel but is excluded from the index because its transcript peaks between the arms (Ukai-Tadenuma et al. 2011); genes differing gene by gene without an index difference stay NO TIME CONFOUND DETECTED with the genes listed and a note that the pattern is not that of a time shift; NOT ASSESSABLE with fewer than 5 usable genes, an arm missing, or a test that cannot reach alpha (≤ 3 units per group; the Welch t on the index and the per-gene table are reported for information). Reason: the plan said "any clock-gene difference → warning" with a panel-level permutation test on the sum of squared t statistics; the true IFN-β claims on Kang showed uniform suppression of clock genes by IFN-β (DBP padj 0.012, CRY1 0.037, PER1 0.051 with 8 vs 8 samples; index p = 0.38), and under that rule they would have been false alarms. The index is biologically motivated (a time shift moves the two arms in opposite directions), but the rule was changed after seeing those results.
- Rhythmic-gene overlap from the Wucher et al. 2023 S8 table (CC BY 4.0), gene-level: per gene the day and night tissue counts and a binomial P, no per-tissue membership. Reason: the plan said tissue-aware; the public table does not allow it, so the check is not tissue-specific and says so.
- Verdict integration: Mode A, imbalance and the effect removed by adjustment → claim NOT SUPPORTED; every other POSSIBLE TIME CONFOUND (shrunk, not estimable, did not hold, Mode B coherent shift) → capped at FRAGILE; NO TIME CONFOUND / NOT ASSESSABLE → unchanged. Documented in HOW_IT_WORKS.
- Bench: the planted time signal goes into all cells of each Kang donor, both IFN conditions, with each donor one sample (`sample = [replicate]`: 8 units, 4 morning vs 4 afternoon); cos-shaped shifts of the panel genes present in the object (11 of 13; CIART and NPAS2 are absent) and of SQSTM1 (night-higher in 6 of 6 GTEx tissues in S8), claimed in Mode B and in Mode A; the true IFN-β claims run with Mode B to measure false alarms. Reason for all cells rather than the control cells of the plan: a donor's ctrl and stim halves were processed together and share one collection time, so pooling them makes the donor the unit and keeps IFN status balanced within every donor, where it cannot confound the fake groups. Reason for the gene: RPS26 was tried first and its donor-to-donor spread in Kang (log2 range about 2.6 across control donors) masked a 0.4 log2 planted shift; SQSTM1's range is about 0.4 within control monocytes. The `TIME_SIGNAL` peak hours are rounded canonical phases (E-box arm morning, PER2 midday, CRY1 delayed to the evening, activator arm evening); they agree with the S8 day/night direction for 11 of 14 genes and disagree for CRY1 (day-higher in 20 tissues), NR1D1 (night-higher in 22) and NR1D2 (no consistent direction); S8 is post-mortem tissue split by time of death, so phases are not directly comparable, and the disagreement is recorded in the code and checked by a test rather than tuned away.

v0.5 status (2026-11-11):
- Done: QC cap with duplicate-point dedupe; stratified within-sample AUROC with `judged_by` / `marker_basis`; candidate rule, per-marker majority, `marker_logic`, Jaccard on kept cells, `reference_is_labelled`; cluster-level doublet diagnostic shown in the report; Scrublet rate fallback; paired tests and `~ pair + condition`; both-groups gate with the not-testable count in the verdict; time check Modes A and B with units, floor, hours rule, estimability and the CRY1-free index; rhythmic-gene table shipped; report sections and methods text; `Thresholds.time_imbalance_hours`.
- Done: benchmark time cases (Mode B rule: label before the check ROBUST, time check POSSIBLE TIME CONFOUND, final verdict not ROBUST; Mode A rule: held before adjustment, effect removed, NOT SUPPORTED) and the strict planted-doublet rule. Measured on the planted object with the current code: balance permutation p = 0.0286 (the floor for 4 vs 4, one unit above where 0.05 becomes reachable) with a +5.8 h circular mean difference; SQSTM1 pseudobulk log2FC 0.66 / padj 0.007 unadjusted and −0.55 / 0.17 after adjustment (direction flips, effect removed); Mode B index permutation p = 0.0286 with 0 of 11 panel genes individually different after BH.
- Done: determinism tests that run the engine twice and compare the full provenance apart from timestamps and runtimes, for a population claim and for DE claims with and without collection times.
- Documented as limitations: a condition-specific population cannot be tested as a composition shift; with ≤ 3 units per group the time permutation tests have no power (Mode A falls back to the hours rule, Mode B is NOT ASSESSABLE); Mode B separates a uniform clock-gene change from a time shift only through the antiphase index; the Scrublet expectation assumes one lane per scored subset.

## v0.6 plan (2026-11-11): populations made of doublets, criteria fixed before measuring

Problem. The planted-doublet case is still a miss. Scrublet flags none of the 150 injected
CD14 monocyte + B cell doublets because they are each other's neighbours; they form their own
cluster only under Pearson residuals (3 of 23 runs), and there the claim holds with doublet
removal on. Goal: a doublet step that removes a population made of heterotypic doublets when
`doublet_removal` is on, without removing real populations.

Candidates, each prototyped and measured before any is adopted:
1. Population-level two-sample test: every cluster of a run against simulated doublets (summed
   raw counts) of candidate parent clusters of different cell types; cross-validated classifier
   AUROC; a cluster that cannot be told apart from its parents' doublets is removed.
2. Cluster-level density of simulated heterotypic doublets: Scrublet-style neighbourhoods with
   simulations drawn only across clusters and library size kept, aggregated per cluster.
3. Per-cell classifier in the style of scDblFinder (Germain et al., F1000Research 2021):
   cluster-based artificial doublets, gradient boosting on neighbourhood, library-size and
   expression features; called cells removed.

Development data: the existing planted-doublet case (PBMC3k + 150 CD14 Mono + B doublets,
seed 0); every cluster of the benchmark grids on Kang 2018 (singlets by demuxlet, 21 runs) and
on unmodified PBMC3k (23 runs); the two Kang GEO lanes with demuxlet calls.

Adoption criteria, fixed now:
- Threshold for a cluster-level statistic: the midpoint between the most doublet-like value
  among Kang and unmodified PBMC3k clusters and the least doublet-like value among planted
  doublet clusters (>= 80 % injected, >= 20 cells), used only if the gap is at least 0.10 on
  the statistic's scale. A per-cell candidate fixes its call rule before it is measured.
- Kang (singlets): over the 21 runs, no cluster removed in at least 19 runs, and the new step
  removes at most 1 % of cells in any run beyond what the current Scrublet rule removes.
- Planted development case: the doublet cluster is removed in every run where it forms with
  `doublet_removal` on.
- GEO lanes: every removed cluster is at least 50 % demuxlet doublets (ambiguous excluded).
- If no candidate meets all four, none is adopted and the attempt is reported as measured.

Held-out evaluation, built only after a candidate and its threshold are fixed, run once and
reported as measured. Two planted cases join the benchmark with the strict rule (caught if the
verdict is not ROBUST, the claim held in at least one run, and every holding run had doublet
removal off):
- H1: 150 doublets of NK + CD14 Mono (seed 11); claim: a GNLY+ NKG7+ CD14+ LYZ+ population exists.
- H2: 150 doublets of CD8 T + B (seed 12); claim: a CD8A+ CD8B+ MS4A1+ CD79A+ population exists.
A held-out block that never forms a population in any run is reported as uninformative (the
strict rule scores it as missed). Non-regression: the 18 existing cases keep their outcome; the
planted-doublet development case may change and is reported separately as a development result.

### v0.6 round 1 (measured 2026-11-11): no candidate adopted

All clusters of the cached grids were measured (Kang 427, PBMC3k 259, planted 250, GEO 115) and
each result was re-computed independently.
- Candidate 1 (two-sample test, parents = single subclusters): least doublet-like planted
  cluster AUROC 0.838 vs most doublet-like real cluster 0.912 (15-cell platelet cluster); gap
  0.074 < 0.10. Seed noise is of the same size: with other seeds the gap was 0.12-0.13 and the
  outcome flipped. Not adopted.
- Candidate 2 (simulated-doublet density): planted clusters tie with a 14-cell NK/T cluster of
  PBMC3k (gap 0). Not adopted. Per cell it ranks doublets better than Scrublet (GEO AUROC 0.88
  vs 0.76 on the baseline run; planted 0.97 vs 0.89).
- Candidate 3 (per-cell classifier, call rule p > 0.5): calls 9 % of Kang singlets beyond
  Scrublet (limit 1 %). Not adopted. Per cell it also ranks better than Scrublet on GEO
  (AUROC 0.91 vs 0.86-0.88; recall 0.82 vs 0.53-0.55 at lower precision).
Diagnostic from the planted data (development data): with whole cell types as parents instead of
single subclusters the planted clusters' AUROC falls to 0.39-0.54, the value a doublet population
should get, because a centroid-chosen subcluster parent (e.g. the high-depth part of the
monocytes) gives simulated doublets larger than the real ones.

### v0.6 round 2 (fixed before measuring): candidate 1 with cell-type groups as parents

- Groups: the other clusters of the run joined by average-linkage clustering of their centroids
  (per-gene z-scored log expression, 1 - Pearson correlation, cut at 1.0), as in round 1.
- Parents: pairs of different groups; the 3 pairs whose group centroids best explain C's centroid
  as a two-part mixture; simulated doublets sum a random cell of each group.
- Statistic: mean over 5 seeded repetitions (new simulations and new folds each time) of the
  out-of-fold AUROC, minimum over the 3 pairs. Everything else as round 1.
- Same threshold rule and criteria; in addition the outcome must be the same for three master
  seeds (0, 1, 2), since round 1 showed seed flips.
- If round 2 fails, nothing is adopted; the per-cell ranking gain of candidates 2 and 3 is
  recorded for later work and the planted-doublet case stays a miss.

### v0.6 round 2 (measured 2026-11-11): adopted

Measured on every cluster of >= 10 cells (Kang 427, PBMC3k 251, planted 245, GEO 115) for master
seeds 0, 1 and 2, and re-implemented independently. Most doublet-like real cluster: AUROC
0.94-0.96 (a 33-cell PBMC3k dendritic-cell cluster); Kang minimum 0.99; GEO minimum 0.99.
Planted doublet clusters: 0.62-0.71. Gap 0.25-0.28 (>= 0.10 required); the per-seed midpoints
are 0.80, 0.83 and 0.82. With any of them: no cluster removed on Kang, GEO or unmodified PBMC3k;
the planted doublet cluster removed in both runs where it forms with doublet removal on.
Implementation threshold, fixed now: 0.80 (the seed-0 midpoint, rounded down, i.e. the most
conservative of the three). The held-out cases H1 and H2 are built and run only after the
implementation is finished.

### v0.6 status (2026-11-11)

- Implemented as fixed in round 2 (threshold 0.80) and committed before the held-out cases ran.
- Benchmark, 20 cases: cases 1-18 keep their outcomes; the step removed no real cluster.
  Development case: the Pearson doublet clusters are removed, but the claim now holds in the
  `n_pcs=20` run (a 293-cell monocyte cluster, 42 % doublets, statistic 0.88) because the
  marker-only claim's reference shifts to that run; still a miss.
- Held-out H1 (NK + CD14 Mono): the block clusters on its own under Pearson residuals with
  statistic 0.95; not removed; missed. H2 (CD8 T + B): merged into a 757-cell real cluster at
  resolution 0.3; missed. Not used to change anything; new held-out cases would be needed.

## Milestone 4 status (2026-11-11): done

- CI on push (`.github/workflows/ci.yml`): ruff check + format, mypy, pytest; offline tests only.
  It was pushed in a separate commit (807ef35) once the token had the `workflow` scope.
- `docs/`: `index.html` (static list of the 16 benchmark reports with verdicts), `reports/`,
  `HOW_IT_WORKS.md`, `report-screenshot.png` (headless Chrome render of a Kang report),
  `.nojekyll`. Built by `scripts/build_docs.py` from `bench/results.json`; committed, no build
  step on Pages. Enable with Settings → Pages → Source "Deploy from a branch", branch `main`,
  folder `/docs`.
- README rewritten with the required sections; Kang quickstart uses `recluster bench --part kang`
  or `scripts/download_kang.py` + `examples/claims/kang_isg15_cd14_monocytes.yaml`.

## Milestone 3 status (2026-11-11): live API call not yet run

- `recluster/llm/`: `common.py` (client, `.env`, number extraction), `parse.py` (claim →
  spec via structured output, data check, one feedback round), `summarize.py` (numeric digest →
  paragraph with a numbers-only guard). CLI: `parse`, `summarize`, `--summarize`.
- 13 tests with fake clients cover the data description, validation feedback, the number guard,
  `.env` loading and the no-key paths. No API key was available in the build environment, so
  the two live calls (`messages.parse` with the `ProposedClaim` schema, `messages.create` for
  the summary) have not been exercised against the real API. First live run should confirm the
  schema is accepted; if the structured-output endpoint rejects a field, the fallback is to
  simplify `ProposedClaim` (all fields are already plain strings, lists and literals).

### Milestone 3 decisions (2026-11-11)

- The model never sees counts, only `.obs` column names and values (up to 20 per column) and
  which gene symbols from the claim text exist in `var_names`. Genes it names that do not exist
  are caught by `check_against` and fed back.
- Structured output via `client.messages.parse(output_format=ProposedClaim)`: a flat schema of
  strings and lists (empty string = not applicable) converted to `ClaimSpec`, rather than
  exposing the nested `ClaimSpec` schema with unions to the API.
- Summaries are generated from a text digest that contains every number the summary may use
  (including derived counts such as runs failed). Numbers glued to letters (gene names) are not
  counted. A violating summary is retried once with the offending numbers named, then discarded.
- Default model `claude-opus-5-5`, the current general-purpose model; adaptive
  thinking is the model default and is left on; no server-side fallback is configured (the
  prompts are low-risk and the layer is optional).
