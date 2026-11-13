# How Recluster works

A plain-language account of every design choice and every statistical test, in the order the
code runs them. File names in parentheses point at the implementation. Default thresholds are
quoted inline; all of them live in `Thresholds` and `GridConfig` (`recluster/config.py`) and are
printed in every report.

## 1. The problems Recluster is built for

A single-cell analysis is a long chain of choices: which cells to discard, how to normalise,
how many genes and components to keep, how finely to cluster, whether to correct for batch,
which random seed to use. Each choice is defensible on its own, and each can change the
answer. A claim that only holds under one particular chain of choices is not a finding about
the biology; it is a finding about the pipeline. Recluster takes a claim and re-runs the
analysis under a bounded set of alternative choices, then reports how often, and under which
choices, the claim survives.

Two failures are documented in the literature and not checked by default in common workflows.
**Pseudoreplication**: cells from one donor are not independent observations of that donor's
condition; treating ten thousand cells as ten thousand replicates makes almost any difference
"significant" (Squair et al., Nat Commun 2021; Zimmerman et al., Nat Commun 2021).
**Time-of-day confounding**: human gene expression changes with the time of collection, and
core clock genes differ between day and night across many tissues (Wucher et al., PLoS Biol
2023); if one condition was sampled in the morning and the other in the afternoon, a "condition
effect" can be a time effect. Recluster's verdicts use the biological replicate as the unit of
inference, as pseudobulk methods do (Crowell et al. 2020), show the naive cell-level result next
to it, and run a time-of-day check on every comparison claim.

The grid is a multiverse or specification-curve analysis (Steegen et al. 2016; Simonsohn et al.
2020) applied to a single claim. Pipeline benchmarks such as pipeComp (Germain et al. 2020) rank
methods; Recluster instead asks whether one claim survives, with replicate-aware tests and a
stated verdict.

## 2. Claims (`recluster/spec.py`)

A claim is a YAML file validated by pydantic. Three types:

- `population_exists`: a population defined by marker genes and/or an existing label exists as a
  distinct cluster.
- `composition_shift`: the population's share of cells differs between two conditions.
- `differential_expression`: named genes change in the population between two conditions.

The spec names the `.obs` columns that matter: `condition`, `sample` (one column or a list that is
combined into one replicate id), `batch` (optional), `labels` (optional) and `collection_time`
(optional, hours or `HH:MM`). An optional `tissue` is recorded with the rhythmic-gene check,
which cannot use it (section 8). `target.marker_logic` is `any` (markers are alternative evidence
for one type; the score is their mean) or `all` (the population co-expresses every marker; the
score is the minimum). `check_against` verifies that every column, value and gene exists before
anything runs, that collection times parse, and that each sample belongs to exactly one
condition. If a donor appears in both conditions (Kang 2018: each patient was split into a
treated and an untreated half) the replicate must be `[donor, condition]`.

### Claim file fields

Keys not listed here are rejected, so a misspelt key stops the run instead of being ignored.
Copy a file from `examples/claims/` (`TEMPLATE.yaml` matches `data/kang_2018.h5ad` as written).

| field | what it is |
|---|---|
| `claim_type` | `population_exists`, `composition_shift` or `differential_expression`; required |
| `name` | free text; required. Also the output folder name: lowercased, non-alphanumerics to '-', first 60 characters |
| `columns.condition` | `.obs` column with the condition; required for composition_shift and differential_expression |
| `columns.sample` | `.obs` column with the biological replicate, or a list of columns joined into one id (e.g. `[donor, condition]`); required for composition_shift and differential_expression |
| `columns.batch` | optional `.obs` column; enables the batch-correction axis |
| `columns.labels` | optional `.obs` column with cell labels; required with `target.label_value` and in audit mode |
| `columns.collection_time` | optional `.obs` column with each sample's collection time: hours, HH:MM, HH:MM:SS or datetime |
| `target.markers` | marker genes of the population; `target` needs `markers` and/or `label_value` |
| `target.negative_markers` | optional genes the population should not express |
| `target.marker_logic` | `any` (default: mean marker score) or `all` (co-expression: minimum score) |
| `target.label_value` | a value of `columns.labels` that names the population |
| `target.min_cells` | smallest matched population that counts, default 20 |
| `comparison.group_a`, `comparison.group_b` | two different values of `columns.condition`; "a" is the group the direction refers to |
| `comparison.direction` | differential_expression: `up_in_a`, `down_in_a` or `either`; composition_shift: `higher_in_a`, `lower_in_a` or `different` |
| `genes` | genes tested; differential_expression only |
| `require` | `all` or `any` of the genes must hold; differential_expression only, default all |
| `tissue` | optional; recorded with the day-night gene overlap, not used (section 8) |

`--grid` and `--thresholds` take YAML files in the same style, with the field names of
`GridConfig` and `GridAxes` (sections 4 and 5) and of `Thresholds` (section 15). Keys left out
keep their defaults, unknown keys are rejected, and every report prints the values used. For
example, three random combinations and two resolutions, and the default alpha written out:

```yaml
# grid.yaml, passed as --grid grid.yaml
n_random: 3
axes: {resolution: [0.5, 1.0]}
```

```yaml
# thresholds.yaml, passed as --thresholds thresholds.yaml
alpha: 0.05
```

## 3. Input modes (`recluster/data.py`)

Recluster always re-derives results from raw counts: `layers['counts']`, then `.raw.X`, then
`.X`, refusing to run if none holds integer counts. `run` mode takes any such object. `audit`
mode is for an object another tool or agent produced: it requires the raw counts stored next to
the processed matrix and a labels column, because an audit checks someone else's labels against
a fresh derivation.

## 4. The pipeline one run executes (`recluster/pipeline.py`)

Every grid point is a `PipelineParams` object with nine fields. One run does, in order:

1. **Cell QC.** Keep cells with at least `min_genes` detected genes and at most `max_mito_pct`
   percent mitochondrial counts. A threshold that alone would remove more than
   `qc_max_cell_loss` (25 %) of the cells is not a defensible alternative on that dataset: it is
   pulled back to the value that removes exactly 25 % (a quantile of the per-cell distribution),
   and the run records both values and labels itself "capped". On Kang 2018, whose median cell
   has 519 detected genes, `min_genes=500` would delete 44 % of cells and becomes 441. Genes
   detected in fewer than 3 cells are dropped. If the matrix has no mitochondrial genes the mito
   axis is removed from the grid rather than silently duplicating runs.
2. **Doublet scoring.** Scrublet (Wolock et al. 2019; scanpy's implementation), run separately
   per sample when every sample has at least 100 cells, with a fixed seed and exact neighbours.
   Scrublet picks its own threshold from the valley of the simulated-doublet histogram, and that
   threshold can land above every observed score, calling nothing (it did on one of the two
   original Kang lanes).
   When the automatic threshold calls fewer than `scrublet_fallback_trigger` (25 %) of the
   expected doublets, the top expected fraction by score is called instead. The expectation is
   `scrublet_expected_rate_per_1000` (0.8 %) per 1,000 cells *of the scored subset*. That equals
   the 10x loading rule only when each sample was captured on its own lane; for samples
   demultiplexed from a pooled lane the true expectation is the lane's. Measured against
   demuxlet on the two Kang GEO lanes with sample = lane, the expectation is 10.8 to 11.1 %
   against 11.5 % called by demuxlet. This stage is cached per QC setting; `doublet_removal`
   then keeps or drops the predicted doublets.
3. **Normalisation and feature selection.** Either total-count scaling to 10,000 per cell,
   `log1p`, selection of `n_hvgs` highly variable genes (Seurat flavour), and z-scaling clipped at
   10; or analytic Pearson residuals (Lause et al. 2021) on the top `n_hvgs` genes.
4. **PCA** with `n_pcs` components (ARPACK solver, seeded).
5. **Batch correction**, optional: Harmony (Korsunsky et al. 2019) on the PCs, keyed by the
   batch column.
6. **Neighbour graph**: 15 nearest neighbours, exact search, so that, with fixed seeds, two runs
   on the same machine give identical provenance apart from timestamps and run times (floats
   rounded to six significant digits).
7. **Leiden clustering** (Traag et al. 2019) at `resolution`, igraph flavour, 2 iterations, seeded.
8. **Doublet populations**, when `doublet_removal` is on and `doublet_cluster_removal` is set (the
   default). Scrublet scores cells one at a time and misses doublets that are each other's
   neighbours, for example many doublets of the same two cell types. This step asks of every
   cluster of at least 10 cells whether it can be told apart from doublets of other cell types.
   The other clusters are joined into cell-type groups (average-linkage on their centroids,
   1 - Pearson correlation of per-gene z-scored log expression, cut at 1.0); the three pairs of
   groups that best explain the cluster's centroid as a two-part mixture are the candidate
   parents; for each, up to 500 doublets are simulated by summing the raw counts of a random cell
   of each group, and a logistic regression on 20 principal components plus library size and
   genes detected is asked to separate the cluster's cells from them (5-fold out-of-fold AUROC,
   averaged over 5 repetitions). The statistic is the lowest AUROC over the three pairs; a cluster
   below `doublet_cluster_auroc_max` (0.80) is removed and the run records which, with its
   statistic and parents. Two held-out cases were defined before the method and threshold were
   fixed, and run once after the step was built. How it was chosen and how it did is in
   `bench/RESULTS.md` ("Doublet populations"): no real cluster removed among 793 measured, the
   development doublet block removed where it forms its own cluster, and both held-out blocks
   missed.

Only the cluster assignment per cell is kept. Expression used by the evaluators is computed once
from the raw counts and is the same in every run.

## 5. The grid (`recluster/grid.py`)

| axis | options | baseline |
|---|---|---|
| min_genes | 100, 200, 500 (capped per dataset, see 4.1) | 200 |
| max_mito_pct | 5, 10, 20 | 10 |
| doublet_removal | off, on | on |
| normalization | log_total, pearson | log_total |
| n_hvgs | 1000, 2000, 3000 | 2000 |
| n_pcs | 20, 30, 50 | 30 |
| resolution | 0.3, 0.5, 1.0, 1.5, 2.0 | 1.0 |
| batch_correction | off, on (only with a batch column) | on |
| seed | 0, 1, 2 | 0 |

The design is the baseline, every one-at-a-time change (so a failure is attributable to one
axis), and 6 seeded random combinations (interactions). A full factorial would be 4,860 runs
without the batch axis and 9,720 with it.

Two grid points whose QC thresholds cap to the same effective values would run the same
pipeline twice and count one configuration twice in the survival rate. Before anything runs,
the engine (`recluster/engine.py`) drops such duplicates, keeping the first of each (the baseline
always), and records every dropped point with the point it duplicates in the provenance
(`merged_runs`, counted in `grid_counts`).

## 6. Finding the population in each run (`recluster/matching.py`)

A claim is about a population; each run has its own clusters. Recluster decides, run by run,
which cells are "the population".

- **Marker score.** For each cell, the mean (or, with `marker_logic: all`, the minimum of the
  z-scores) of CP10k-log1p expression of the marker genes, minus the mean of the negative
  markers, z-scored across all cells. Computed once from raw counts.
- **Marker enrichment, judged within samples.** For a cluster, a stratified within-sample AUROC
  of the marker score: every sample that holds at least one cell of the cluster and at least one
  other cell contributes the Mann-Whitney U of the cluster's cells against the other cells *of
  the same sample*, and the AUROC is the sum of those U statistics divided by the sum of
  n1 × n0 (cluster cells times other cells) over the contributing samples. A sample's weight is
  therefore its number of cluster-versus-other pairs; there is no per-sample size requirement.
  A cluster of at least `min_cluster_size` (10) cells qualifies at AUROC `marker_auroc_min`
  (0.80). Cells are only ever ranked against cells of their own sample, which is what matters: a
  QC setting that removes lymphocytes preferentially makes the "rest" of the dataset
  myeloid-rich, and an AUROC against all other cells then drops for a real monocyte cluster
  (IFN-stimulated monocytes keep only one of three classic markers). The within-sample value is
  immune to that, and it also rejects dendritic-cell clusters that the whole-dataset AUROC
  admitted through LYZ alone. The statistic was first implemented as the minimum over samples
  and abandoned after measuring it on simulations: a population with true AUROC 0.85 and 12
  cells in each of 16 samples was rejected 95 % of the time by the minimum and 0 % by the
  stratified statistic (mean 0.851), and at true AUROC 0.90 the minimum still rejected about 32 %,
  because one sample's few cells happening to rank low was enough. A cluster with no contributing
  sample, and every cluster when there is no sample column, is judged by the whole-dataset
  AUROC. Each run records which statistic judged each candidate (`judged_by`) and sums it up
  for the run as `marker_basis`: within-sample, overall or mixed.
- **Per-marker majority.** With two or more markers, more than `marker_majority_min` (50 %) of
  them must be individually enriched in the cluster (per-gene AUROC at least
  `marker_per_gene_auroc_min`, 0.60). A marker set carried by a single gene, or a
  "monocyte-and-B-cell" set satisfied by monocytes alone, does not qualify. With
  `marker_logic: all`, every marker must be enriched.
- **Label overlap.** When a label is given, a cluster qualifies when at least
  `label_precision_min` (50 %) of its cells carry it.
- **Candidate rule.** When a claim gives both markers and a label, a cluster qualifies by either
  criterion (`candidate_rule: either`; `markers` or `labels` restricts it). Qualifying by labels
  needs a reference built from the claim's label column. A reference carried over from the
  baseline run's match (section "Reference set" below) is the thing under test, so it is used
  only for the Jaccard overlap and never as a candidate criterion: otherwise a merged
  low-resolution cluster holding most of the baseline population would qualify "by labels" with
  no marker evidence at all. On Kang, CD4 T cells are split by condition into two clusters whose
  marker evidence falls short in some runs; the label criterion recovers them. The matched
  population is the union of qualifying clusters, and its whole-dataset AUROC against the other
  kept cells is reported for every run. For population claims that union AUROC is checked again
  (section 7.1); comparison claims apply no union-AUROC check beyond candidacy.
- **Reference set and Jaccard.** The labelled cells when a label is given; otherwise the matched
  set of the first run, normally the baseline, in which the matched population reaches marker
  AUROC 0.80 (recorded as `reference_source`). Jaccard is computed over barcodes that passed the
  run's QC, and the fraction of reference cells that QC removed is reported separately, so "QC
  deleted part of it" and "it is not there" are distinguished.
- **Sample spread.** Samples contributing at least 5 cells, and the largest share held by one
  sample.
- **Doublet diagnostics (informational).** The matched population's median Scrublet score next
  to the dataset's, and a mixture test asking whether its centroid lies on the segment between
  two other clusters with a larger library size (`evaluators/doublet_cluster.py`). Both are
  shown in the report's population panel and never change a verdict: no cutoff tried on the
  benchmark's real clusters and its planted doublet clusters separated true mixtures from
  genuine intermediate populations with a margin worth enforcing, so the two thresholds
  (`doublet_residual_ratio_max`, `doublet_libsize_ratio_min`) only label the diagnostic.

## 7. Evaluators (`recluster/evaluators/`)

### 7.1 population_exists

A run holds when the matched set has at least `min_cells` (20), its marker AUROC (the
whole-dataset AUROC of the matched set against the other kept cells) is at least 0.80 (if
markers were given), its Jaccard to the reference is at least `jaccard_min` (0.50), and it has at
least 5 cells in each of at least `min_samples_present` (2) samples.

### 7.2 composition_shift

For each sample: `n` cells in the run, `k` of them in the population, and the logit of a
smoothed proportion, `log((k + 0.5) / (n - k + 0.5))` (the logit transformation is one of the
two that propeller offers; Phipson et al. 2022). When the same donors appear in both
groups (`evaluators/pairing.py` detects this from the replicate key: every pair level must have
a sample in each group) the decisive test is a paired t-test on the per-donor differences of
those logits, with a Wilcoxon signed-rank test alongside; otherwise a Welch t-test on the logits
with a Mann-Whitney U test alongside. Pairing is the correct model for a split-sample design and
removes donor baseline variance. Measured on Kang with this code, over its 8 labelled cell types
(7 testable): the unpaired Welch test calls none (all p ≥ 0.67); the paired test calls one,
CD14+ monocytes (p = 0.035; about 23.8 % of cells in stim against 25.4 % in ctrl, lower in stim
for 7 of 8 donors), and standard errors shrink 2 to 9 fold. That cell type is the benchmark's
null composition case, and the result is reported as measured.

A run is testable only when the population has at least 5 cells in at least 2 samples of *each*
group. When the matching finds the population in one condition only, the tool cannot tell a
biological absence from the clustering having lost that condition's version of the population,
and the first benchmark's only composition false alarm arose exactly this way: a population
matched in one condition produced a 0.1 % versus 35 % "shift". The cost of the rule is that a
genuinely condition-specific population cannot be tested as a composition shift; such runs are
reported as not testable, they count as failures, and the verdict states how many runs were not
testable. The claim holds when p < 0.05 and the direction matches. A chi-square on the pooled
cells is shown for contrast only.

### 7.3 differential_expression

Two tests per gene on the cells of the matched population.

- **Cell-level**: two-sided Wilcoxon rank-sum on CP10k-log1p values, log2 fold change of means
  with a pseudocount of 1, Benjamini-Hochberg over the claimed genes.
- **Replicate-aware (pseudobulk)**: raw counts summed per sample within the population (samples
  with fewer than 5 cells dropped; genes with fewer than 10 total counts dropped, claimed genes
  kept); DESeq2 (Love et al. 2014, through PyDESeq2; Muzellec et al. 2023) with design
  `~ condition`, or `~ pair + condition` when every pair has a pseudobulk in both groups, Wald
  test. P-values of the claimed genes are BH-adjusted over
  the claimed genes; the genome-wide adjusted value is shown as well. A run is testable only with
  at least 2 pseudobulk replicates per group; otherwise it is reported as not testable and counts
  as a failure, and the verdict says how many runs were.
- A gene holds at either level when the adjusted p is below 0.05 and the fold change has the
  claimed direction; the claim holds when all (or, with `require: any`, any) of its genes hold at
  the replicate level. Each gene is tagged `both`, `cell_only`, `pseudobulk_only` or `neither`.
- **Pseudoreplication flag.** Raised when the cell-level test supports the claim in at least 80 %
  of runs while the replicate-aware test supports it in fewer than 50 %.

## 8. The time-of-day check (`recluster/time_confound.py`)

Runs once per comparison claim, after the grid, on the baseline run's matched population (or on
the labelled cells when the baseline found nothing). Never at the level of the cell.

**Units.** The unit of both permutation tests is the sample, unless a column of the sample key
other than the condition column groups several samples per level, every level sits in one
condition (and, when times are given, at one collection time), and each group holds at least 2
levels; the coarsest such column is the unit. Example: a condition assigned per donor with
`sample = [donor, label]` makes the donor the unit, because a donor's two samples are not two
independent observations of its condition or of when it was sampled. The report names the unit
column when one was used.

**The permutation tests.** A label-symmetric statistic is compared with its value under every
reassignment of the group labels to the units when there are at most 20,000 assignments,
otherwise under 5,000 seeded draws. The smallest p the exact test can produce, 2/C(n, k) for
equal groups and 1/C(n, k) otherwise, is reported with every result. With three or fewer units
per group no permutation test can reach 0.05 (3 vs 3 gives 20 assignments and a floor of 0.10;
4 vs 4 gives 70 and a floor of 0.0286, one unit above the point where 0.05 becomes reachable).

**Mode A, collection time known** (`columns.collection_time`). Each sample's time is the circular
mean of its cells' times (a spread above half an hour within a sample is reported), and each
unit's time becomes an angle.

*Balance.* The statistic is the distance between the groups' mean (cos, sin) vectors; the circular
mean of each group and their difference in hours are reported, with Welch t-tests on cos and sin
alongside. The groups count as imbalanced when the permutation p is below `alpha` (0.05) **or**
the circular means differ by at least `time_imbalance_hours` (2.0 h), and the report says which
rule applied. The hours rule exists because of the floor: with 3 vs 3 units the test cannot go
below 0.10, and a 2 h shift moves a rhythmic gene by up to 2·sin(π·2/24) ≈ 0.52 of its
half-amplitude, 26 % of its peak-to-trough swing, which should not pass as "balanced".

*Adjustment* (only when imbalanced). The claim's replicate-level model is refitted with the same
design the verdict used (`~ condition`, or `~ pair + condition` when every pair has a pseudobulk
in both groups) plus cos(2πt/24) and sin(2πt/24), with DESeq2's Cook's outlier filtering off so
every claimed gene keeps a p-value in both fits; composition claims use an ordinary least-squares
fit of the logit proportions on condition, with pair fixed effects when the design is paired,
with and without the two covariates. A time column that adds no rank is dropped first: with two
distinct clock times sin is an affine function of cos, and one time per pair is absorbed by the
pairing. When the condition coefficient of the adjusted model is not estimable (too few
replicates for the terms, no time term left, or collection time collinear with the condition,
for example every control sample at one hour and every treated sample at another) the report
says the effect cannot be separated from time, and the verdict is capped at FRAGILE, never NOT
SUPPORTED. Otherwise, for each claimed gene, the report shows the unadjusted and adjusted
estimate and adjusted p (BH over the claimed genes in each fit), the shrinkage, whether the
direction flips and whether the gene survives (adjusted p < 0.05, same sign, claimed direction),
combined over genes with the claim's `require`.

*Decision.* Balanced → NO TIME CONFOUND DETECTED. Imbalanced and the claim did not hold anyway →
POSSIBLE TIME CONFOUND (the adjustment is uninformative; the verdict is capped). Imbalanced and
not estimable → POSSIBLE (capped). Imbalanced, the claim held, and after adjustment it is no
longer significant **and** the estimate shrinks by at least `time_shrink_flag` (50 %) or flips →
POSSIBLE with the effect removed → NOT SUPPORTED. No longer significant but the estimate barely
moved → POSSIBLE (capped: the adjusted model spends degrees of freedom on time, and time and
condition cannot be separated with this design). Still significant but shrunk by at least 50 % →
POSSIBLE (capped). Otherwise → NO TIME CONFOUND DETECTED. "Held before adjustment" is the
baseline run's own outcome from the grid, so the check never contradicts the grid.

**Mode B, collection time unknown.** Per-unit log2 CPM, from all cells of the unit, of a panel
of 13 canonical core clock genes: ARNTL (BMAL1), CLOCK, NPAS2, PER1, PER2, PER3, CRY1, CRY2,
NR1D1, NR1D2, DBP, TEF, CIART. This is the transcriptional-translational feedback loop as
reviewed by Takahashi (Nat Rev Genet 2017;18:164-179), and the genes Wucher et al. (2023) report
as day-night variable across human tissues. A gene is usable with at least
`time_min_gene_counts` (10) counts across units, and the mode is assessable with at least
`time_min_panel_genes` (5) usable genes, both arms of the clock represented in the index, and
at least 2 units per group. Two statistics are computed. Per-gene Welch t-tests across units
with BH correction answer the literal question "do clock genes differ between the groups". The
decisive statistic is a **clock-state index** per unit: the mean z-score (across units) of the
E-box arm, PER1, PER2, PER3, CRY2, NR1D1, NR1D2, DBP, TEF and CIART, minus the mean z-score of
the activator arm, ARNTL, CLOCK and NPAS2, compared between groups with the permutation test
on the absolute mean difference. CRY1 stays in the per-gene panel but is left out of the index
because its transcript peaks several hours after the other E-box targets, between the two arms
(Ukai-Tadenuma et al., Cell 2011), and would blur an index built on antiphase. The index is
decisive because a collection-time difference moves the two arms in opposite directions, whereas
a treatment that suppresses clock genes moves them together. When the index differs
(permutation p < 0.05) the label is POSSIBLE TIME CONFOUND with a warning that a
collection-time or circadian confound is possible; when genes differ gene by gene but the index
does not, the label stays NO TIME CONFOUND DETECTED, the report lists the genes and says the
pattern is not that of a time shift. When the test cannot reach 0.05 (three or fewer units per
group) Mode B is NOT ASSESSABLE; the reason names the units and the floor, and the Welch t on
the index and the per-gene table are reported for information.
Mode B is a red-flag screen and does not estimate circadian phase; the report says so.

**How the Mode B rule was arrived at.** The original design was "any clock-gene difference →
warning", with a panel-level permutation test on the sum of squared t statistics. It was
replaced by the index after running the true IFN-β claims on Kang: IFN-β suppresses clock genes
uniformly in the stimulated samples (DBP padj 0.012, CRY1 0.037, PER1 0.051 with 8 vs 8 samples),
and the index does not call it (p = 0.38). Under the original rule those true claims would have
been false alarms. The index is biologically motivated, since a time shift moves the two arms in
opposite directions, but the rule change was made after seeing those results, and the benchmark's
false-alarm count for Mode B should be read with that in mind.

**Rhythmic-gene overlap.** The public S8 Table of Wucher et al. (CC BY 4.0), rebuilt as
`recluster/resources/daynight_genes.csv` by `scripts/fetch_daynight_genes.py`, lists for every
day-night gene the number of GTEx tissues in which it is higher by day or by night. It is
gene-level, so the check is not tissue-specific and a `tissue` given in the spec is only noted;
claimed genes that are day-night in at least 5 tissues with a consistent direction (binomial
P < 0.05) are flagged with their tissue counts.

**Verdict integration** (`recluster/verdict.py`). POSSIBLE TIME CONFOUND in Mode A with the
effect removed: the claim becomes NOT SUPPORTED. Every other POSSIBLE TIME CONFOUND (the effect
shrinks but survives, the design cannot separate time from condition, the claim did not hold
before adjustment, or Mode B finds a phase-coherent clock-gene difference): a ROBUST claim is
capped at FRAGILE. NO TIME CONFOUND DETECTED and NOT ASSESSABLE leave the verdict unchanged; the
report records the verdict before and after.

## 9. Verdicts (`recluster/verdict.py`)

Survival is the fraction of grid runs in which the claim held. ROBUST at 80 % or more; FRAGILE at
20 % or more; NOT SUPPORTED below; then the time-of-day rule above. Runs that could not test the
claim (section 7.2 and 7.3) count as failures, and the verdict summary states how many there
were. The report lists, for every setting, how many runs with that setting failed, and which
one-at-a-time changes flipped the baseline result.

## 10. Determinism and provenance (`recluster/provenance.py`)

Every stochastic step is seeded; neighbour search is exact; permutation tests enumerate or use a
fixed seed. Floating-point results are rounded to six significant digits so two runs compare
equal. `provenance.json` records package versions, platform, the input's SHA-256 and shape, the
spec, the grid, the thresholds, the dropped duplicate grid points, every run's effective QC
values and outputs, and the time check. Tests run the engine twice and compare the full
provenance apart from timestamps and runtimes, for a population claim and for differential
expression claims with and without collection times.

## 11. Reports (`recluster/report/`)

One self-contained HTML file per claim: inline CSS and SVG, no scripts, no external assets. The
verdict is shown with text and an icon, never by colour alone. The population panel shows the
matched population's marker AUROC, Jaccard, sample spread and the two doublet diagnostics; the
replicate-aware panel labels its tests as paired or unpaired according to the design the run
used; the time-of-day section shows the balance test with its floor, the refit, or the clock
panel. The methods paragraph is produced from a fixed template filled with the settings actually
used, including the design the baseline run used. Nothing in the report comes from a language
model unless the optional LLM summary was requested, in which case it is boxed and labelled.

### How to read a report

A report goes from the answer to the evidence. From the top:

- **Verdict.** ROBUST when the claim held in at least 80 % of runs, FRAGILE at 20 % to 80 %, NOT
  SUPPORTED below 20 %, then the time-of-day rule (section 8). The sentence next to the badge
  says how many runs held, whether the baseline held, and, for comparison claims, the
  time-of-day result. "Survival" is
  the share of runs in which the claim held. "Naive cell-level" is how often the cell-level test
  would have called the claim; it is shown for contrast and never decides the verdict.
- **Runs.** A run is one complete re-analysis of the same raw counts. The grid is about 21
  re-analyses: run 0 is the baseline with the default settings, the one-at-a-time runs each
  change one setting, and the random runs change several (section 5). Runs are numbered from 0
  in every table.
- **Pseudoreplication flag**, when shown: the cell-level test supports the claim in at least
  80 % of runs and the donor-level test in fewer than 50 % (section 7.3).
- **Survival across the grid.** One square per run; the thick-bordered square is the baseline.
  "By analysis choice" shows, for each value of each setting, the share of runs with that value
  in which the claim held; ● marks the baseline value.
- **What broke it.** First the one-at-a-time changes that flipped the baseline result; read
  these first, because each differs from the baseline in one setting only. Then "Failure rate by
  setting": for each value, how many runs used it and how many failed. When the baseline fails,
  the baseline values (●) show high failure rates simply because they occur in most runs, so they
  do not single out a cause.
- **Replicate-aware vs cell-level.** For comparison claims, the baseline run's numbers for both
  tests (per gene, or per sample for composition). For a population claim, the matched
  population in the baseline run and whether it is present in several samples rather than
  driven by one.
- **Time-of-day check**, for comparison claims: the balance test with its floor and the refit
  when collection times are given, or the clock-gene panel when they are not (section 8).
- **As labelled**, when the claim gives a label: the same test on the cells exactly as labelled,
  without re-clustering; this is the claim as the original analysis made it.
- **All runs**, **Methods**, **Thresholds** and **Provenance**: every run's settings and
  outcome, the methods paragraph, the cutoffs used, and the inputs and package versions.

What a verdict does not mean, and what to do next:

- **ROBUST** means the claim survives the grid of analysis choices, tested across donors; it does
  not mean the claim is true. A flaw built into the experimental design survives every analysis
  choice: in Kang 2018 each condition was captured on its own 10x lane, so a lane artifact would
  also come out ROBUST (section 14). Next: check the design, and read "What broke it" for
  settings that failed together, such as the seeds.
- **FRAGILE** means the claim depends on the analysis choices, or that a possible time-of-day
  confound capped the verdict. It does not mean the claim is false. Next: read "What broke it"
  for the settings it depends on, and the time-of-day section; report the claim with those
  settings.
- **NOT SUPPORTED** means the claim held in fewer than 20 % of runs, or that the effect
  disappeared once collection time was adjusted for. It does not show that the opposite is
  true. Next: check the pseudoreplication flag, and in "All runs" whether the population was
  found at all; runs that could not test the claim count as failures, and the verdict summary
  says how many there were.

## 12. What the benchmark taught (`bench/RESULTS.md`)

The first benchmark run (v0.1) showed the grid's absolute QC values were not adapted to the
data, that whole-dataset marker AUROCs move with the cell-type mix, that an additive marker score
cannot express co-expression, and that Scrublet's automatic threshold can fail outright. Each
fix in sections 4 to 7 came from one of those misses, was measured on the data before it was
adopted, and the benchmark was re-run with both tables kept.

The planted-doublet case is scored with a strict rule: it counts as caught only if the verdict is
not ROBUST, the claim held in at least one run, and every run in which it held had doublet
removal off. A claim that never holds says nothing about doublet removal and is not counted as
caught. The cluster-level mixture check is shown in the report's population panel as a
diagnostic and takes no part in that rule.

The time-of-day benchmark (`recluster/bench/kang.py`) plants a signal into *all* cells of each
Kang donor, both IFN conditions included, so IFN status is balanced within every donor and cannot
confound the fake groups; each donor is one sample (`sample = [replicate]`), giving 8 units,
4 vs 4. The 8 donors are split by a seeded shuffle into a morning group (08:00, 09:00, 10:30,
12:00) and an afternoon group (13:00, 15:00, 16:30, 18:00). For the 11 panel clock genes present in the object (CIART and NPAS2 are absent) and one output
gene the counts of every cell of a donor are multiplied by 2^(a·cos(2π(t − peak)/24)) and
Poisson-resampled, so the injected effect is a smooth function of time rather than a group
label. The output gene is SQSTM1, higher by night in 6 of 6 GTEx tissues in the S8 table (0 day,
6 night). Which gene carries the signal was a design decision, not a tuning of the signal: RPS26
(day-higher in 16 of 16 tissues) was tried first and rejected because its donor-to-donor spread
in Kang, a log2 range of about 2.6 across the control donors, masked a 0.4 log2 planted shift;
SQSTM1 has a range of about 0.4 within control monocytes. The peak hours in `TIME_SIGNAL` are
rounded canonical phases (E-box arm in the morning, PER2 around midday, CRY1 delayed to the
evening, activator arm in the evening) with half-amplitudes of 0.2 to 0.7 log2. Checked against
the S8 day/night direction they agree for 11 of 14 genes and disagree for CRY1 (day-higher in 20
tissues in S8), NR1D1 (night-higher in 22) and NR1D2 (no consistent direction); the disagreement
is recorded in the code and checked by a test rather than tuned away. The S8 table comes from
post-mortem GTEx tissue classified by time of death, so its day/night split and the planted peak
hours are not directly comparable.

Measured with the current code on the planted object: the Mode A balance test gives permutation
p = 0.0286, the smallest value 4 vs 4 can produce, with a +5.8 h circular mean difference
(afternoon minus morning). SQSTM1 has a pseudobulk log2 fold change of 0.66 and padj 0.007
before adjustment; after adjusting for cos/sin of collection time the fold change is −0.55 with
padj 0.17, so the direction flips and the effect is removed, and the verdict is NOT SUPPORTED. In
Mode B the index permutation p is 0.0286 while 0 of the 11 usable panel genes differ individually
after BH: the index sees what single genes do not. The true IFN-β claims run with Mode B to
measure false alarms (section 8).

## 13. The LLM layer and its boundary (`recluster/llm/`)

`recluster parse` turns a sentence into a claim spec: the model sees the `.obs` columns and
values and which gene symbols from the sentence exist in the data, returns a structured object,
and the result is validated and checked against the file before the user confirms it.
`--summarize` writes a paragraph from a numeric digest of a finished run; a guard rejects any
paragraph containing a number that is not in the digest. The model never sees counts and never
produces a number that appears in a verdict, a table or a chart.

## 14. Known limitations

- ROBUST means a claim survives the grid of analysis choices, tested across donors; it does not
  mean the claim is true. A flaw built into the experimental design survives every analysis
  choice: in Kang 2018 each condition was captured on its own 10x lane, so a lane artifact would
  also come out ROBUST. Apart from collection time, Recluster cannot detect confounds built into
  the design.
- The survival fraction depends on how many grid runs vary each setting. On the 21-run Kang
  grid three runs use a seed other than the baseline's, so a claim that fails in all three still
  holds in 18 of 21 runs and is ROBUST. The report's table of what broke the claim lists such
  settings; read it with the verdict.
- The Mode A adjustment fits the same 24-hour cosine that the planted benchmark signal follows.
  The benchmark shows the mechanism works when its model is right, not how it performs on real
  collection-time effects, which need not follow a single 24-hour cosine.
- The clock-gene screen pools all cells of each unit into one pseudobulk and does not adjust for
  cell-type composition. It has been tested on a planted signal and on three group comparisons
  on real data without a planted time effect (stim vs ctrl and two fake donor groupings), not on
  real time-stamped single-cell data.
- Apart from the planted SQSTM1 shift (log2 fold change 0.66), the effects the benchmark
  recovers are large (ISG15 in CD14+ monocytes, log2 fold change 7.71 with the labels as given);
  sensitivity to modest real effects is not measured.
- Each false-by-construction benchmark case uses one seeded donor split; the false-positive rate
  across splits is not estimated.
- The within-sample AUROC needs a sample column; without one, matching falls back to the
  whole-dataset AUROC with its sensitivity to cell-type mix.
- A genuinely condition-specific population cannot be tested as a composition shift: the
  both-groups rule reports such runs as not testable, because the tool cannot tell a biological
  absence from a clustering that lost one condition's version of the population.
- With three or fewer units per group the time permutation tests have no power: Mode A falls back
  to the hours rule, and Mode B is NOT ASSESSABLE.
- Mode B cannot distinguish a treatment that changes clock genes uniformly from a collection-time
  difference except through the antiphase index; a treatment that happened to move the two arms
  in opposite directions would be flagged as a possible time confound.
- The public day-night list is gene-level, so the overlap check is not tissue-specific.
- The Scrublet expectation assumes one lane per scored subset; for samples demultiplexed from a
  pooled lane the fallback calls fewer doublets than the lane holds.
- Composition tests assume one sample per donor and condition; nested designs beyond that are
  not modelled.
- The mixture-based doublet diagnostic has no reliable cutoff and is reported, not enforced.
- Doublet detection is the weakest step: Scrublet flags almost none of a block of same-type
  doublets, and the doublet-population step missed both held-out blocks: the NK + monocyte block
  formed its own cluster and was not recognised (statistic 0.95), and the CD8 T + B block never
  formed a cluster of its own. Doublets merged into a real cluster are beyond any cluster-level
  test.

## 15. Glossary of thresholds (defaults)

| name | default | used by |
|---|---|---|
| robust_min_survival | 0.80 | verdict |
| fragile_min_survival | 0.20 | verdict |
| alpha | 0.05 | all tests |
| marker_auroc_min | 0.80 | candidate clusters (stratified within-sample AUROC), population recovery (whole-dataset AUROC of the matched set) |
| marker_per_gene_auroc_min | 0.60 | per-marker majority |
| marker_majority_min | 0.50 | per-marker majority (strictly more than) |
| candidate_rule | either | markers and/or labels (labels only with a reference from the label column) |
| jaccard_min | 0.50 | population recovery |
| jaccard_on_kept_cells | true | Jaccard denominator |
| label_precision_min | 0.50 | candidate clusters in label mode |
| min_cluster_size | 10 | candidate clusters |
| min_samples_present | 2 | replicate presence, both-groups rule |
| min_cells_per_sample | 5 | sample spread, pseudobulk inclusion |
| pseudobulk_min_total_counts | 10 | genes entering DESeq2 |
| composition_test | logit_t | decisive composition test |
| doublet_residual_ratio_max | 0.60 | cluster-level doublet diagnostic (informational, labels only) |
| doublet_libsize_ratio_min | 1.50 | cluster-level doublet diagnostic (informational, labels only) |
| time_min_gene_counts | 10 | usable clock-panel gene |
| time_min_panel_genes | 5 | Mode B assessable |
| time_shrink_flag | 0.50 | Mode A shrinkage flag |
| time_imbalance_hours | 2.0 | Mode A: circular mean difference that counts as imbalanced |
| qc_max_cell_loss (grid) | 0.25 | QC threshold cap |
| scrublet_expected_rate_per_1000 (grid) | 0.008 | Scrublet fallback (per scored subset) |
| scrublet_fallback_trigger (grid) | 0.25 | Scrublet fallback |
| doublet_cluster_removal (grid) | true | doublet-population step on when doublet_removal is on |
| doublet_cluster_auroc_max (grid) | 0.80 | doublet-population step: remove clusters below this AUROC |
| pseudoreplication_cell_survival | 0.80 | flag |
| pseudoreplication_replicate_survival | 0.50 | flag |

## 16. References

- Crowell, H. L. et al. muscat detects subpopulation-specific state transitions from
  multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11, 6077 (2020).
  https://doi.org/10.1038/s41467-020-19894-4
- Germain, P.-L., Lun, A., Macnair, W. & Robinson, M. D. Doublet identification in single-cell
  sequencing data using scDblFinder. *F1000Research* 10, 979 (2021).
  https://doi.org/10.12688/f1000research.73600.1
- Germain, P.-L., Sonrel, A. & Robinson, M. D. pipeComp, a general framework for the evaluation
  of computational pipelines, reveals performant single cell RNA-seq preprocessing tools.
  *Genome Biol* 21, 227 (2020). https://doi.org/10.1186/s13059-020-02136-7
- Kang, H. M. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic
  variation. *Nat Biotechnol* 36, 89-94 (2018). https://doi.org/10.1038/nbt.4042
- Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with
  Harmony. *Nat Methods* 16, 1289-1296 (2019). https://doi.org/10.1038/s41592-019-0619-0
- Lause, J., Berens, P. & Kobak, D. Analytic Pearson residuals for normalization of single-cell
  RNA-seq UMI data. *Genome Biol* 22, 258 (2021). https://doi.org/10.1186/s13059-021-02451-7
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2. *Genome Biol* 15, 550 (2014).
  https://doi.org/10.1186/s13059-014-0550-8
- Muzellec, B., Teleńczuk, M., Cabeli, V. & Andreux, M. PyDESeq2: a python package for bulk
  RNA-seq differential expression analysis. *Bioinformatics* 39, btad547 (2023).
  https://doi.org/10.1093/bioinformatics/btad547
- Phipson, B. et al. propeller: testing for differences in cell type proportions in single cell
  data. *Bioinformatics* 38, 4720-4726 (2022). https://doi.org/10.1093/bioinformatics/btac582
- Simonsohn, U., Simmons, J. P. & Nelson, L. D. Specification curve analysis. *Nat Hum Behav* 4,
  1208-1214 (2020). https://doi.org/10.1038/s41562-020-0912-z
- Squair, J. W. et al. Confronting false discoveries in single-cell differential expression.
  *Nat Commun* 12, 5692 (2021). https://doi.org/10.1038/s41467-021-25960-2
- Steegen, S., Tuerlinckx, F., Gelman, A. & Vanpaemel, W. Increasing transparency through a
  multiverse analysis. *Perspect Psychol Sci* 11, 702-712 (2016).
  https://doi.org/10.1177/1745691616658637
- Takahashi, J. S. Transcriptional architecture of the mammalian circadian clock. *Nat Rev
  Genet* 18, 164-179 (2017). https://doi.org/10.1038/nrg.2016.150
- Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected
  communities. *Sci Rep* 9, 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z
- Ukai-Tadenuma, M. et al. Delay in feedback repression by cryptochrome 1 is required for
  circadian clock function. *Cell* 144, 268-281 (2011). https://doi.org/10.1016/j.cell.2010.12.019
- Wolock, S. L., Lopez, R. & Klein, A. M. Scrublet: computational identification of cell
  doublets in single-cell transcriptomic data. *Cell Syst* 8, 281-291.e9 (2019).
  https://doi.org/10.1016/j.cels.2018.11.005
- Wucher, V., Sodaei, R., Amador, R., Irimia, M. & Guigó, R. Day-night and seasonal variation of
  human gene expression across tissues. *PLoS Biol* 21, e3001986 (2023).
  https://doi.org/10.1371/journal.pbio.3001986
- Zimmerman, K. D., Espeland, M. A. & Langefeld, C. D. A practical solution to pseudoreplication
  bias in single-cell studies. *Nat Commun* 12, 738 (2021).
  https://doi.org/10.1038/s41467-021-21038-1
