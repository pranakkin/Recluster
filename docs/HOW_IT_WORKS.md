# How Recluster works

A plain-language account of every design choice and every statistical test, in the order the
code runs them. File names in parentheses point at the implementation. Default thresholds are
quoted inline; all of them live in `Thresholds` (`recluster/config.py`) and are printed in every
report.

## 1. The problem Recluster is built for

A single-cell analysis is a long chain of choices: which cells to discard, how to normalise,
how many genes and components to keep, how finely to cluster, whether to correct for batch,
which random seed to use. Each choice is defensible on its own, and each can change the
answer. A claim that only holds under one particular chain of choices is not a finding about
the biology; it is a finding about the pipeline. Recluster takes a claim and re-runs the
analysis under a bounded set of alternative choices, then reports how often, and under which
choices, the claim survives.

The second problem it is built for is pseudoreplication. Cells from one donor are not
independent observations of that donor's condition. Treating ten thousand cells as ten thousand
replicates makes almost any difference "significant". Recluster's verdicts use the biological
replicate (donor, sample, mouse) as the unit of inference and show the naive cell-level result
next to it so the gap is visible.

## 2. Claims (`recluster/spec.py`)

A claim is a YAML file validated by pydantic. Three types:

- `population_exists`: a population defined by marker genes and/or an existing label exists as a
  distinct cluster.
- `composition_shift`: the population's share of cells differs between two conditions.
- `differential_expression`: named genes change in the population between two conditions.

The spec names the `.obs` columns that matter: `condition`, `sample` (one column or a list that is
combined into one replicate id), `batch` (optional) and `labels` (optional). `check_against`
verifies that every column, value and gene exists before anything runs, and that each sample
belongs to exactly one condition. That last check matters: if a donor appears in both conditions
(as in Kang 2018, where each patient was split into treated and untreated halves) the replicate
must be `[donor, condition]`, otherwise the test would compare a donor with itself.

Why YAML and pydantic: the spec is the contract between the human (or the LLM) and the engine.
It must be readable, diffable and rejected loudly when wrong.

## 3. Input modes (`recluster/data.py`)

Recluster always re-derives results from raw counts. It looks for integer counts in
`layers['counts']`, then `.raw.X`, then `.X`, and refuses to run if none is found. `run` mode
takes any such object. `audit` mode is for an object another tool or agent produced: it requires
the raw counts to be stored alongside the processed matrix and requires a labels column, because
the point of an audit is to check someone else's labels against a fresh derivation.

## 4. The pipeline one run executes (`recluster/pipeline.py`)

Every grid point is a `PipelineParams` object with nine fields. One run does, in order:

1. **Cell QC.** Keep cells with at least `min_genes` detected genes and at most `max_mito_pct`
   percent mitochondrial counts (gene names starting with `MT-` or `mt-`). Genes detected in
   fewer than 3 cells are dropped. If the matrix has no mitochondrial genes (Kang 2018 does not)
   the mito axis is removed from the grid rather than silently running duplicate points.
2. **Doublet scoring.** Scrublet (scanpy's implementation), run separately per sample when every
   sample has at least 100 cells, because doublets form within a 10x lane. The seed is fixed and
   neighbours are exact, so the score is a function of the QC setting only; this stage is cached
   per QC setting. When `doublet_removal` is on, predicted doublets are removed.
3. **Normalisation and feature selection.** Either total-count scaling to 10,000 per cell,
   `log1p`, selection of `n_hvgs` highly variable genes (Seurat flavour), and z-scaling clipped at
   10; or analytic Pearson residuals on the top `n_hvgs` genes. These are the two mainstream
   choices and they disagree most for lowly expressed genes.
4. **PCA** with `n_pcs` components (ARPACK solver, seeded).
5. **Batch correction**, optional: Harmony on the PCs, keyed by the batch column, when one is
   declared and the axis is on.
6. **Neighbour graph**: 15 nearest neighbours, exact search (sklearn). Approximate search is
   faster but not bitwise reproducible across machines, and reproducibility is a design goal.
7. **Leiden clustering** at `resolution`, igraph flavour, 2 iterations, seeded.

Only the cluster assignment per cell is kept. Expression used by the evaluators is computed
once from the raw counts and is the same in every run, so differences between runs come only
from which cells were kept and how they were grouped.

## 5. The grid (`recluster/grid.py`)

The default axes and their options:

| axis | options | baseline |
|---|---|---|
| min_genes | 100, 200, 500 | 200 |
| max_mito_pct | 5, 10, 20 | 10 |
| doublet_removal | off, on | on |
| normalization | log_total, pearson | log_total |
| n_hvgs | 1000, 2000, 3000 | 2000 |
| n_pcs | 20, 30, 50 | 30 |
| resolution | 0.3, 0.5, 1.0, 1.5, 2.0 | 1.0 |
| batch_correction | off, on (only with a batch column) | on |
| seed | 0, 1, 2 | 0 |

The design is the baseline, plus every one-at-a-time change (each non-baseline value of each
axis, 16 or 17 runs), plus 6 seeded random combinations drawn from the full product. One-at-a-time
runs make failures attributable: if the only change was `resolution=0.3` and the claim fell, that
is the cause. Random combinations probe interactions. A full factorial would be 4,860 runs without the batch
axis and 9,720 with it; the bounded design keeps a small dataset under a few minutes and a
25,000-cell dataset at roughly five to eight minutes per claim.

The values are meant to be defensible, not extreme; whether they are is itself a finding. The
benchmark showed that `min_genes=500` is not a defensible alternative on a dataset whose median
cell has 519 genes (section 11).

## 6. Finding the population in each run (`recluster/matching.py`)

A claim is about a population; each run has its own clusters. Recluster has to decide, run by
run, which cells are "the population".

- **Marker score.** For each cell, the mean of CP10k-log1p expression of the marker genes minus
  the mean of the negative markers, z-scored across all cells. Computed once from raw counts.
- **Marker enrichment per cluster.** The AUROC of the marker score for cells in a cluster versus
  all other cells, computed through the Mann-Whitney U identity. AUROC is scale-free and does
  not depend on cluster size, unlike a mean score.
- **Candidate clusters (marker mode).** Clusters with AUROC at least 0.80 and at least 10 cells.
  Their union is the matched set. Taking a union tolerates a population that a high resolution
  split into pieces; a population merged into a larger cluster lowers the AUROC and fails, which
  is the intended meaning of "exists as a distinct population".
- **Candidate clusters (label mode).** When only a label is given: clusters in which at least
  50 % of cells carry the label, again at least 10 cells.
- **Reference set.** The labelled cells when a label is given; otherwise the matched set of the
  baseline run (or of the first run that finds it). The report says which.
- **Jaccard** between the matched set and the reference, over barcodes, counting cells a run
  removed by QC as absent. This is deliberate: a QC setting that deletes half the population
  should lower the score. It is also what made the Kang NK claims FRAGILE under `min_genes=500`
  (section 11), and the honest reading is that both are true.
- **Sample spread.** The number of samples contributing at least 5 cells, and the largest share
  of the population held by one sample. A "population" that comes from one donor is a batch
  effect until proven otherwise.

Marker mode is used whenever markers are given, even if a label is also given; the label then
defines the reference. The reason is that the point of re-derivation is to find the population
without trusting the labels.

## 7. Evaluators (`recluster/evaluators/`)

### 7.1 population_exists

A run holds when the matched set has at least `min_cells` (20), its marker AUROC is at least
0.80 (if markers were given), its Jaccard to the reference is at least 0.50, and it has at least
5 cells in each of at least 2 samples (when a sample column exists). The replicate condition is
deliberately weak; it is a sanity check, not a test.

### 7.2 composition_shift

For each sample: `n` cells in the run, `k` of them in the population, and the logit of a
smoothed proportion, `log((k + 0.5) / (n - k + 0.5))`. The pseudo-count keeps empty and full
samples finite. The decisive test is a Welch t-test on those logits between the samples of
group A and the samples of group B (unequal variances; 8 vs 8 donors in Kang). A Mann-Whitney U
test on the raw proportions is reported next to it and can be made decisive in the thresholds
file. The effect is the difference in mean logit (a log odds ratio) and the ratio of mean
proportions. A run holds when p < 0.05 and the direction matches the claim.

The naive alternative, a chi-square test on the pooled 2 × 2 table of cells, is computed only
for contrast. It treats every cell as a replicate; on Kang it gives p = 0.008 for a difference
the replicate-aware test puts at p = 0.67.

With fewer than two replicates per group the test is undefined, and the run fails with a warning
rather than falling back to a cell-level test.

### 7.3 differential_expression

Two tests per gene, on the cells of the matched population split by condition.

- **Cell-level**: a two-sided Wilcoxon rank-sum (Mann-Whitney U) test on CP10k-log1p values, and
  a log2 fold change of mean normalised expression with a pseudocount of 1. P-values are
  Benjamini-Hochberg adjusted over the claimed genes.
- **Replicate-aware (pseudobulk)**: raw counts are summed per sample within the population;
  samples with fewer than 5 cells are dropped; genes with fewer than 10 total counts are dropped
  (claimed genes are always kept); DESeq2 (pydeseq2, design `~ condition`, Wald test) estimates a
  log2 fold change and p-value per gene. The p-values of the claimed genes are BH-adjusted over
  the claimed genes, because a claim names its genes in advance; the genome-wide adjusted value
  is shown as well.
- A gene holds at either level when the adjusted p is below 0.05 and the fold change has the
  claimed direction. The claim holds when all (or, with `require: any`, any) of its genes hold
  at the replicate level. Each gene is tagged `both`, `cell_only`, `pseudobulk_only` or `neither`.
- **Pseudoreplication flag.** Raised for composition and expression claims when the cell-level
  test supports the claim in at least 80 % of runs while the replicate-aware test supports it in
  fewer than 50 %. The reference is Squair et al., Nature Communications 12:5692 (2021), which
  showed that cell-level tests on multi-sample data produce false discoveries at high rates and
  that pseudobulk methods do not.

Why pseudobulk rather than a mixed model: with a handful of donors, a mixed model is fragile and
slow; summing counts per donor and using a count model built for a handful of replicates is the
better-understood option, and it is what the Squair benchmark found to perform best.

## 8. Verdicts (`recluster/verdict.py`)

Survival is the fraction of grid runs in which the claim held. ROBUST at 80 % or more; FRAGILE at
20 % or more; NOT SUPPORTED below. The report also lists, for every setting, how many runs with
that setting failed, and which one-at-a-time changes flipped the baseline result. The thresholds
are round numbers chosen before any benchmark ran and have not been changed since.

## 9. Determinism and provenance (`recluster/provenance.py`)

Every stochastic step is seeded: Scrublet, PCA, Harmony, Leiden, and the random grid draws. The
neighbour search is exact. Floating-point results in the provenance log are rounded to six
significant digits so two runs can be compared with equality. `provenance.json` records package
versions, platform, the input file's SHA-256 and shape, the spec, the grid, the thresholds, and
for every run its parameters, cell and cluster counts, runtime and evaluator outputs. A test
runs the engine twice and asserts identical provenance apart from timestamps.

## 10. Reports (`recluster/report/`)

One self-contained HTML file per claim: inline CSS and inline SVG, no scripts, no external
assets, so it can be emailed or served from a static folder. The verdict is shown with text and
an icon, never by colour alone. The methods paragraph is produced from a fixed template filled
with the settings actually used; the sentence about mitochondrial QC changes when the data has
no mitochondrial genes. Nothing in the report comes from a language model unless the optional
LLM summary was requested, in which case it is boxed and labelled.

## 11. What the benchmark taught (`bench/RESULTS.md`)

The full table and every miss are in `bench/RESULTS.md`. The short version:

- The replicate-aware tests behave: all four false claims built from real data are NOT
  SUPPORTED, while the cell-level tests would have called three of them in 81 to 100 % of runs.
- The grid's absolute QC values are not adapted to the data. On Kang 2018, a shallow dataset,
  `min_genes=500` removes 44 % of cells and over half of the NK cells, and every true claim came
  out FRAGILE because of that one setting. The report attributes it correctly, but a relative
  threshold (a quantile of genes per cell) would be a better alternative to offer.
- Marker-driven matching has edges. IFN-β changes monocyte marker expression enough that the
  stimulated monocyte cluster's AUROC lands at 0.79 under tighter QC, just under the 0.80 cutoff.
  When that happens for a composition claim, the population "disappears" from one condition and
  the test calls a shift that is not there. A "present in both groups" precondition is the fix.
- An additive marker score cannot express co-expression, so a claim about a doublet population
  (monocyte and B-cell markers together) matched real monocytes and B cells. A co-expression
  score is the fix. Scrublet also missed a block of 150 identical synthetic doublets, because
  they are each other's nearest neighbours; and against demuxlet on the real Kang lanes its
  score ranks doublets well (AUROC 0.87) while its automatic threshold called none in one lane.

None of these were patched before publishing the numbers; they are listed as proposed changes.

## 12. The LLM layer and its boundary (`recluster/llm/`)

Two optional conveniences, both labelled as LLM output. `recluster parse` turns a sentence into
a claim spec: the model sees the `.obs` columns and values and which gene symbols from the
sentence exist in the data, returns a structured object, and the result is validated and
checked against the file before the user confirms it. `--summarize` writes a paragraph from a
numeric digest of a finished run; a guard rejects any paragraph containing a number that is
not in the digest. The model never sees counts and never produces a number that appears in a
verdict, a table or a chart.

## 13. Known limitations

- Existence is penalised for QC removal (section 6) and the QC axes are absolute (section 5).
- Matching is global across conditions; a treatment that remodels a population can make it
  drop out of one condition (section 11).
- Doublet removal depends on Scrublet's automatic threshold, which can fail.
- Multi-marker claims use an additive score.
- The replicate check for existence claims is a presence count, not a test.
- Runs within a claim are sequential; only the benchmark parallelises across claims.
- Composition tests are unpaired even when the design is paired (Kang's donors appear in both
  conditions). A paired test would have more power; it was not implemented in v0.

## 14. Glossary of thresholds (defaults)

| name | default | used by |
|---|---|---|
| robust_min_survival | 0.80 | verdict |
| fragile_min_survival | 0.20 | verdict |
| alpha | 0.05 | composition and expression tests |
| marker_auroc_min | 0.80 | candidate clusters, population recovery |
| jaccard_min | 0.50 | population recovery |
| label_precision_min | 0.50 | candidate clusters in label mode |
| min_cluster_size | 10 | candidate clusters |
| min_samples_present | 2 | replicate check for existence |
| min_cells_per_sample | 5 | sample spread, pseudobulk inclusion |
| pseudobulk_min_total_counts | 10 | genes entering DESeq2 |
| composition_test | logit_t | decisive composition test |
| pseudoreplication_cell_survival | 0.80 | flag |
| pseudoreplication_replicate_survival | 0.50 | flag |
