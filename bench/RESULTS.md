# Benchmark results

Measured on 2026-11-11 with Recluster 0.5.0, default grid and default thresholds, on an Apple
M5 Pro laptop (`recluster bench --jobs 4`, 25 minutes wall time). Every number below is copied
from `results.json`; the per-run numbers come from the `provenance.json` next to each report
under `bench/out/`. Where the outcome is a miss or a false alarm, it says so.

This is the second run of the benchmark. The first run (v0.1, 16 cases) is kept as
`results_v0.1_table.md` and `results_v0.1.json`, and its narrative is summarised in the last
section; the changes made between the two runs, each with its reason and the measurement that
motivated it, are listed under "What changed between the runs". The pass rule of every case
was written before the case ran. One rule was tightened between the runs (the planted-doublet
case, see note C) and two new cases were added (the planted time-of-day effect, note D).

Verdict rule for all cases: ROBUST if the claim held in >= 80 % of grid runs, FRAGILE if in
>= 20 %, otherwise NOT SUPPORTED, then the time-of-day rule (a possible confound caps ROBUST at
FRAGILE; an effect that the time adjustment removes makes the claim NOT SUPPORTED). "Held"
counts runs where the replicate-aware test passed. "Cell-level" is the share of runs where the
naive cell-level statistic (Wilcoxon on cells, or chi-square on pooled cells) would have called
the claim; it never enters the verdict. "Time check" is the label of the time-of-day check:
`none` = NO TIME CONFOUND DETECTED, `POSSIBLE` = POSSIBLE TIME CONFOUND; population claims
have no comparison and no time check.

## Headline table

| # | Case | Expected | Verdict | Held | Cell-level | Time check | Outcome |
|---|---|---|---|---|---|---|---|
| 1 | Kang: ISG15 up in CD14+ monocytes, stim vs ctrl | ROBUST | ROBUST | 21/21 | 100 % | none | correct |
| 2 | Kang: ISG15, IFI6, IFIT3, MX1 all up in CD14+ monocytes | ROBUST | ROBUST | 21/21 | 100 % | none | correct |
| 3 | Kang: ISG15 up in CD4 T cells | ROBUST | ROBUST | 21/21 | 100 % | none | correct |
| 4 | Kang: CD14+ monocytes exist across donors | ROBUST | ROBUST | 21/21 | – | – | correct |
| 5 | Kang: NK cells exist across donors | ROBUST | ROBUST | 20/21 | – | – | correct (note A) |
| 6 | Kang: CD14+ monocyte proportion differs stim vs ctrl (no known shift) | NOT SUPPORTED | FRAGILE | 11/21 | 71 % | none | **false alarm** (note B) |
| 7 | False: ISG15 differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 1/21 | 100 % ⚑ | none | correct |
| 8 | False: LYZ differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 100 % ⚑ | none | correct |
| 9 | False: monocyte proportion differs between fake donor groups | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 24 % | none | correct |
| 10 | False: ISG15 up in fake group A (condition shuffled at donor level) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 100 % ⚑ | none | correct |
| 11 | Planted time effect claimed as a condition effect, collection times not given | held before the check, then flagged | FRAGILE (ROBUST before the check) | 21/21 | 100 % | POSSIBLE | caught (note D) |
| 12 | The same planted effect with collection times given | held, then removed by adjustment | NOT SUPPORTED (ROBUST before the check) | 21/21 | 100 % | POSSIBLE | caught (note D) |
| 13 | Planted: 150 synthetic doublets as a fake population | not ROBUST; holds only with doublet removal off | NOT SUPPORTED | 3/23 | – | – | **missed** (note C) |
| 14 | Planted: over-clustering shard of CD4 T cells | not ROBUST | NOT SUPPORTED | 0/23 | – | – | caught |
| 15 | Control: a shard that is a whole population (B cells) | ROBUST | ROBUST | 23/23 | – | – | control OK |
| 16 | Planted: QC threshold deletes NK cells | failure table names max_mito_pct | NOT SUPPORTED | 4/23 | – | – | caught |
| 17 | Planted: batch artifact claimed as a CD4 T subtype | not ROBUST; holds only with batch correction off | NOT SUPPORTED | 4/24 | – | – | caught |
| 18 | Control: the real CD4 T population under the batch effect | ROBUST | ROBUST | 24/24 | – | – | control OK |

⚑ = pseudoreplication flag raised by the report (cell-level statistics call it in >= 80 % of
runs while replicate-aware statistics call it in < 50 %).

Score by the pre-registered rules: 16 of 18 as expected (5 true claims ROBUST, 4 false claims
NOT SUPPORTED, both planted time cases caught, 3 planted errors caught, 2 controls OK); 1 false
alarm (case 6) and 1 miss (case 13). The time-of-day check raised no false alarm: on the 8
comparison claims without a planted time effect its label was NO TIME CONFOUND DETECTED every
time (clock-state index permutation p = 0.38 on the real IFN-β comparisons, 0.63 on the fake
donor splits).

### The same cases in the first run

| # | Case | v0.1 verdict (held) | v0.2 verdict (held) |
|---|---|---|---|
| 1 | ISG15 up in CD14+ monocytes | FRAGILE (16/21) | ROBUST (21/21) |
| 2 | Four ISGs up in CD14+ monocytes | FRAGILE (16/21) | ROBUST (21/21) |
| 3 | ISG15 up in CD4 T cells | FRAGILE (15/21) | ROBUST (21/21) |
| 4 | CD14+ monocytes exist | FRAGILE (16/21) | ROBUST (21/21) |
| 5 | NK cells exist | FRAGILE (15/21) | ROBUST (20/21) |
| 6 | Monocyte proportion, stim vs ctrl (null) | FRAGILE (8/21), false alarm | FRAGILE (11/21), false alarm |
| 7-10 | Four false claims | NOT SUPPORTED (0-1/21) | NOT SUPPORTED (0-1/21) |
| 13 | Planted doublets | ROBUST (23/23), missed | NOT SUPPORTED (3/23), missed |
| 14-18 | Over-clustering, QC deletion, batch effect and their controls | all caught / controls OK | all caught / controls OK |

The five true claims went from FRAGILE to ROBUST because the one grid value that broke them,
`min_genes=500` on a dataset whose median cell has 519 genes, is now capped at the value that
removes a quarter of the cells (441), and because a population is now judged within each sample
rather than against a cell-type mix that the QC setting changes (details below). The null case
is a false alarm in both runs but for a different reason each time, and the doublet case is a
miss in both runs but now with the fake population rejected rather than accepted.

## Notes on the misses and the false alarm

**A. NK cells (case 5) fail in one run of 21.** The run is the random combination
`min_genes=500 (capped to 441), normalization=pearson, n_hvgs=1000, n_pcs=20,
resolution=0.3`: at resolution 0.3 with Pearson residuals and 1,000 features the NK cells sit
in a cluster shared with CD8 T cells, the best cluster's stratified marker AUROC is 0.78
(cutoff 0.80), and the label criterion does not rescue it either (the cluster is less than half
NK). The other 20 runs match 1,700 to 1,800 cells with Jaccard 0.85 to 0.88 against the
`cell_type` labels. In the first run this case failed six runs, five of them from the
`min_genes=500` cell loss alone.

**B. The null composition claim (case 6) is a false alarm by the pre-registered rule, and the
reason is not the one of the first run.** In the first run the 8 holding runs were an artifact:
the marker matching found the monocytes in one condition only, so the per-sample test compared
0.1 % against 35 %. That cannot happen now (the matched population is both clusters in every
run, Jaccard 0.85 to 0.97 against the labels, present in all 16 samples). What holds now is a
small, consistent difference that the paired test sees and the unpaired test did not: with the
labels as given, the CD14+ monocyte share is 25.4 % in the control halves and 23.8 % in the
stimulated halves of the same eight donors, lower in stim for 7 of the 8 (paired t on logit
shares p = 0.035, Wilcoxon signed-rank p = 0.016; the first run's unpaired Welch test gave
p = 0.67). Across the grid the paired p ranges from 0.0002 (`resolution=0.3`) to 0.64, and 11
of 21 runs fall under 0.05, so the verdict is FRAGILE (52 %), with Pearson normalisation and
`min_genes=100` among the settings under which it does not hold. The case was registered as
"no known shift" because the same PBMCs were split and one half stimulated for six hours, which
is a judgment about the biology, not a planted fact; a consistent 1.6-percentage-point drop in
the recovered monocyte share after IFN-β is plausible (activation changes survival, capture and
clustering of monocytes). It is reported as a false alarm because that is what the rule says,
and the rule was not changed after the result. Two things are worth reading off it: the paired
test is the right test for this design and is more sensitive than the unpaired one, and a
FRAGILE verdict with the breakers listed is the tool doing its job on a claim whose truth is not
known.

**C. The planted doublet population (case 13) is missed by the strict rule, while the fake
population is now rejected.** The claim holds in 3 of 23 runs, all three with Pearson-residual
normalisation, under which the 150 injected monocyte+B doublets form their own cluster of 74 to
100 cells with all four markers enriched (per-marker AUROC 0.70 to 0.87, cluster AUROC 0.94 to
0.95). Under log-normalisation they do not form a cluster at all: at the baseline 113 of them
land in the CD14 monocyte cluster and 33 in the B cluster, no cluster carries more than half of
the four markers, and the matched population is empty. Two of the three holding runs had doublet
removal on, so the rule "holds only when doublet removal is off" is not met and the case counts
as a miss: Scrublet (automatic threshold 0.26; fallback not triggered because the automatic
threshold called 1.0 % of cells against an expected 2.3 %, above the 25 % trigger) flagged 0 of
the 150 injected cells, exactly as in the first run. The verdict NOT SUPPORTED is right, but for
the wrong reason: the fake population fails because it exists only under one normalisation, not
because the doublets were removed. The cluster-level mixture diagnostic on the three Pearson
clusters gives residual ratios of 0.58 to 0.68 and library-size ratios of 1.40 to 1.63; the
cutoff analysis below had predicted that no cutoff separates such clusters with margin, and the
report shows the diagnostic without acting on it.

**D. The planted time-of-day effect (cases 11 and 12) is caught in both modes, and the naive
analysis falls for it.** The construction (details below) gives 8 donors, 4 collected between
08:00 and 12:00 and 4 between 13:00 and 18:00, with the 11 core clock genes present in the object and SQSTM1 shifted
in every cell as a cosine of the hour. The claim "SQSTM1 up in CD14+ monocytes in the afternoon
group" holds in 21 of 21 runs: the pseudobulk DESeq2 test on 4 vs 4 donors gives log2 fold
change 0.66 with adjusted p = 0.007, and the cell-level test p = 4e-10. A standard analysis
would report it as a condition effect. With the collection times hidden (case 11) the clock-gene
screen flags it: the clock-state index (E-box arm minus activator arm) separates the two groups
completely, permutation p = 0.029, which is the smallest value 4 vs 4 units can give, while 0 of
the 11 usable panel genes differ individually after BH correction; the label is POSSIBLE TIME
CONFOUND and the verdict is capped from ROBUST to FRAGILE. With the times given (case 12) the
balance test says the groups were collected 5.8 hours apart (permutation p = 0.029, again the
floor) and the refit with cos/sin of collection time turns the fold change from +0.66 (padj
0.007) into −0.55 (padj 0.17): the direction flips, the effect is removed, and the verdict is
NOT SUPPORTED with the reason "the effect is explained by collection time". The day-night gene
table flags SQSTM1 as night-higher in 6 of 6 GTEx tissues where it varies. Both pass rules
require the claim to have held before the check (so that a flag on a claim nobody would have
believed does not count) and both are met. Read with the caveat that 4 vs 4 is one unit above
the point where a permutation test can reach 0.05 at all; with 3 vs 3 donors Mode B would have
been NOT ASSESSABLE and Mode A would have relied on the 2-hour rule.

On the real IFN-β comparisons the screen does not fire. The 8 vs 8 samples differ in two clock
genes gene by gene (DBP, adjusted p 0.012, and CRY1, 0.037, both lower after IFN-β; PER1 at
0.051), but the index does not move (p = 0.38): IFN-β lowers both arms together, which is not
the antiphase pattern of a time shift. The report says so and the verdicts stay ROBUST. How the
Mode B rule was arrived at, and that it was settled after seeing this result, is recorded in
`docs/HOW_IT_WORKS.md` section 8 and in `PLAN.md`.

## What changed between the runs

Each change was measured on the data before it was adopted; the measurement is given with the
reason. The numbers below, except the doublet mixture analysis (item 9), are recomputed by
`scripts/design_measurements.py` (about one minute after `recluster bench --part kang`), which
writes `bench/design_measurements.md` with each cited value next to the measured one.
Thresholds introduced here were fixed before the second run started, and none was
changed after seeing the second run's results.

1. **QC thresholds capped at 25 % cell loss** (`GridConfig.qc_max_cell_loss`). The first run's
   five FRAGILE true claims all failed in the `min_genes=500` runs, which remove 44 % of Kang's
   cells (median 519 genes per cell). The cap pulls the threshold back to the 25 % quantile,
   441 genes, and the run is labelled "capped". Grid points that become identical after capping
   are dropped and recorded; none were in this benchmark.
2. **Marker AUROC judged within samples, stratified.** The first run compared a cluster against
   all other cells, and that comparison moves with the cell-type mix a QC setting keeps: the
   stimulated monocyte cluster scored 0.83 at baseline and 0.79 under `min_genes=500` although
   the cluster itself was the same, because the filter preferentially removes lymphocytes and
   leaves a LYZ-richer background. The statistic is now the Mann-Whitney U summed over samples,
   divided by the summed n1·n0 (every sample with at least one cluster cell and one other cell
   contributes). A first version took the minimum over samples instead; on simulations (16
   samples × 12 cluster cells, true AUROC 0.85) that minimum rejected the population 95 % of
   the time against 0 % for the stratified statistic (mean 0.851), so the minimum was not used.
3. **Per-marker majority.** A cluster qualifies by markers only when more than half of the
   markers are individually enriched (per-gene AUROC >= 0.60). This is what rejects the
   monocyte-or-B reading of the planted doublet claim under log-normalisation. `marker_logic:
   all` (minimum of per-marker z-scores) is available for co-expression claims.
4. **Labels or markers as candidates** (`candidate_rule: either`) when a claim gives a label
   column; a reference carried over from the baseline run is used only for the Jaccard overlap.
   Measured on Kang before adoption: the CD4 T population, which marker enrichment alone loses
   at low resolution (its two condition-split halves each score 0.75 to 0.80 against the rest),
   is recovered in every run with Jaccard 0.90 to 0.97.
5. **Jaccard over the cells a run kept**, with the fraction of reference cells removed by QC
   reported separately, so "QC deleted part of it" is no longer read as "it is not there".
6. **Paired replicate tests** when every donor appears in both groups (paired t on per-donor
   logit differences; `~ pair + condition` in DESeq2). Measured on Kang with the labels as given:
   standard errors shrink 2- to 9-fold, the unpaired Welch test calls none of the 7 testable
   cell types, the paired test calls one, CD14+ monocytes (note B); for DE in monocytes, 95 % of
   the unpaired design's calls are kept under blocking with no sign flips.
7. **A population must be present in at least 2 samples of both groups** for a composition or
   expression test; otherwise the run is not testable, counts as a failure, and the verdict says
   how many runs were. This is what removes the first run's composition artifact (note B of
   that run). The cost, that a genuinely condition-specific population cannot be tested as a
   composition shift, is documented as a limitation.
8. **Scrublet fallback.** When Scrublet's automatic threshold calls fewer than a quarter of the
   expected doublets (0.8 % per 1,000 cells of the scored subset), the top expected fraction by
   score is called. Against demuxlet on the two Kang lanes (table below) the expected rate is
   within 0.7 percentage points of demuxlet's and the fallback recovers the stimulated lane,
   where the automatic threshold had called nothing. It does not catch the planted PBMC3k
   doublets: their scores are low (median 0.10) and even an oracle cutoff at the true count would
   catch 28 % of them.
9. **Cluster-level doublet mixture check: informational only.** The plan had a check that
   explains a population's centroid as a mixture of two other clusters and flags it when the
   fit is good and the library size is high. Measured over 472 real clusters from 32 runs on
   Kang and PBMC3k plus the planted doublet clusters: neither the residual ratio nor the
   library-size ratio separates on its own (15 real clusters fall below the planted clusters'
   worst ratio, 53 exceed their smallest library-size ratio); a conjunction of the two does
   separate, but the feasible residual cutoff lies in (0.647, 0.691), a margin of 0.02 on each
   side. No rule with that margin is defensible, so the diagnostic is shown in the report's
   population panel and never changes a verdict.
10. **The planted-doublet pass rule** is the strict one: caught only if the verdict is not
    ROBUST, the claim held in at least one run, and every run in which it held had doublet
    removal off. A looser reading ("a claim that never holds counts as caught") was written
    while building and withdrawn, because a tool that finds nothing would pass it.
11. **The time-of-day check** (new; `docs/HOW_IT_WORKS.md` section 8). Decisions recorded with
    their reasons: the unit of the permutation tests is the donor when samples are nested in
    donors of one condition (two samples of one donor are not two observations of when it was
    sampled); the smallest attainable permutation p is reported and, because three or fewer
    units per group cannot reach 0.05, the groups also count as imbalanced when their mean
    collection times differ by at least 2 hours (`time_imbalance_hours`); the adjustment refits
    the same design the verdict used, and a design in which time cannot be separated from the
    condition caps the verdict at FRAGILE instead of removing the effect; CRY1 is kept in the
    per-gene panel but left out of the clock-state index because its transcript peaks between
    the two arms; and the index, rather than any single-gene difference, decides Mode B. That
    last rule replaced the plan's "any clock-gene difference warns" after the real IFN-β claims
    showed uniform suppression of clock genes (note D); it is biologically motivated but was
    settled after seeing that result.
12. **The planted time signal.** The claimed gene was changed from RPS26 to SQSTM1 before the
    second run: RPS26's donor-to-donor spread in Kang (a log2 range of about 2.6 across the
    control donors) masked a 0.4 log2 planted shift, while SQSTM1 varies by about 0.4 log2
    within control monocytes. The planted half-amplitude of SQSTM1 was raised from 0.4 to 0.6
    log2 (1.5-fold at the peak; a 2.3-fold peak-to-trough swing, within the range reported for
    rhythmic human genes) after a smoke test showed that 0.4 gave an unadjusted pseudobulk
    padj of 0.065, that is, a signal the naive analysis would not have called. This is a choice
    of the planted signal's strength, made so that requirement (a) of the case, that a standard
    analysis calls the effect, is met; it is recorded here as such. The clock-gene amplitudes
    (0.2 to 0.7 log2) and peak hours were not changed.

## Doublets vs demuxlet (GEO GSE96583, batch 2)

The processed object holds singlets only, so Scrublet was compared with demuxlet on the two
original lanes (29,065 barcodes; 24,366 singlets, 3,169 doublets, 1,217 ambiguous; 313
barcodes without a call). Scrublet ran exactly as in the pipeline: per lane, baseline QC, exact
neighbours, automatic threshold with the rate-based fallback. Ambiguous barcodes are excluded.

| Lane | Cells | demuxlet doublet rate | Called | Precision | Recall | AUROC of score |
|---|---|---|---|---|---|---|
| both | 27,406 | 11.5 % | 10.3 % | 0.60 | 0.54 | 0.87 |
| ctrl | 13,853 | 11.5 % | 9.1 % | 0.67 | 0.53 | 0.88 |
| stim | 13,553 | 11.5 % | 11.5 % (fallback) | 0.55 | 0.55 | 0.86 |

In the first run (automatic threshold only) the stimulated lane had 0.0 % called and recall 0,
and the combined recall was 0.27. demuxlet can only see doublets from two different individuals
(about 7/8 of doublets with 8 donors), so precision is a lower bound. The fallback's expected
rate assumes each scored subset was captured on its own lane, which is true here and not in
general.

## How each case was built

### Real data (Kang et al. 2018, Nat Biotechnol 36:89)

Processed object from pertpy (`https://exampledata.scverse.org/pertpy/kang_2018.h5ad`,
SHA-256 `e6a5adac…9830`): 24,673 cells, 15,706 genes, raw counts; `label` (ctrl 12,315 / stim
12,358), `replicate` (8 donors), `cell_type` (8 types). The object has no mitochondrial genes,
so the `max_mito_pct` axis is dropped and the grid has 21 runs (14 one-at-a-time, 6 random).
No batch column is declared because condition is confounded with 10x lane.

- Cases 1 to 5: unmodified data; replicate = donor × condition (16 samples, paired on the
  donor). Populations are located by markers (CD14/LYZ/S100A8; IL7R/CD3D/CCR7; GNLY/NKG7/KLRD1
  minus CD3D) or by the `cell_type` labels, with the labels as the reference set.
- Case 6: unmodified data; a composition claim with no planted or known effect.
- Cases 7 to 9: control cells only; donors shuffled with seed 0 and split into
  A = patient_1016, patient_1039, patient_107, patient_1256 and
  B = patient_101, patient_1015, patient_1244, patient_1488. Nothing differs between A and B
  by construction.
- Case 10: all cells; each of the 16 donor × condition samples assigned to A or B so that each
  fake group holds 4 ctrl and 4 stim samples (seed 0).
- Cases 11 and 12: all cells of all 8 donors. Donors shuffled with seed 0; the first four
  (patient_1016, patient_107, patient_1039, patient_1256) form group A with collection times
  08:00, 09:00, 10:30 and 12:00, the other four (patient_1244, patient_101, patient_1015,
  patient_1488) form group B with 13:00, 15:00, 16:30 and 18:00. Each donor is one sample: its
  control and stimulated cells are pooled, so IFN-β status is balanced within every donor and
  cannot confound the groups, and the analysis has 8 units (4 vs 4). For each of the 11 core clock
  genes present in the object (CIART and NPAS2 are not) and for SQSTM1, the counts of every cell of a donor
  are multiplied by 2^(a · cos(2π(t − peak)/24)) with the donor's time t and the gene's peak
  hour and half-amplitude a from `TIME_SIGNAL` in `recluster/bench/kang.py`, then Poisson
  resampled. The injected effect is therefore a smooth function of the hour, not of the group
  label. SQSTM1 peaks at 20:00 with half-amplitude 0.6 log2 and is higher by night in 6 of 6
  GTEx tissues where it varies (Wucher et al. 2023, S8). The peak hours of the clock genes are
  rounded canonical phases (E-box arm in the morning, PER2 around midday, CRY1 delayed to the
  evening, activator arm in the evening); checked against the S8 day/night direction they agree
  for 11 of 14 genes and disagree for CRY1, NR1D1 and NR1D2, a disagreement that is recorded in
  the code and checked by a test rather than tuned away (S8 is post-mortem tissue classified by
  time of death, so its day/night split and the planted phases are not directly comparable).
  Case 11 passes `sample: [replicate]` and `condition: fake_group` without times; case 12 adds
  `collection_time`.

### Planted errors (PBMC3k, `scanpy.datasets.pbmc3k()`, 2,700 cells)

Reference labels come from Recluster's baseline pipeline at resolution 1.0, named by markers
(CD4 T 1,164; CD14 Mono 450; B 332; CD8 T 313; FCGR3A Mono 186; NK 160; DC 31; Mk 15;
49 filtered). Grid: 23 runs (24 when a batch column exists). Seed 0 throughout.

- Case 13: 150 synthetic doublets, each the summed raw counts of a random CD14 monocyte and a
  random B cell, appended as new cells. Claim: a population with markers CD14, LYZ, MS4A1,
  CD79A exists. Caught if not ROBUST, the claim held in at least one run, and every holding run
  had doublet removal off.
- Cases 14 and 15: Leiden at resolution 3.0 gives 30 clusters; CD4 T cells are split into 16
  shards. The largest shard (cluster_14, 116 cells, 100 % CD4 T) is claimed as a population via
  labels in audit mode. Control: cluster_2 (328 cells, Jaccard 0.99 with B cells). Caught if
  not ROBUST; control passes if ROBUST.
- Case 16: mitochondrial counts of the 160 reference NK cells rescaled so each cell's mito
  fraction is uniform in 12 to 16 % (all other cells stay below 6 %). The baseline
  `max_mito_pct=10` removes all of them. Claim: NK cells (GNLY, NKG7, KLRD1; not CD3D) exist.
  Caught if the top row of the failure table is `max_mito_pct`. Result: the claim held only in
  the 4 runs with `max_mito_pct=20`; the top breaker is `max_mito_pct`.
- Cases 17 and 18: a random 40 % of cells assigned to "batch2"; 300 expressed genes multiplied
  by 2.5 and 300 others divided by 2.5 in those cells. Clustering without correction yields
  a 450-cell cluster (98 % batch2, CD4 T), labelled "CD4 T subtype" and claimed in audit mode
  with `batch` declared. Caught if not ROBUST and every holding run had batch correction off
  (the 4 holding runs all had it off). Control: the reference CD4 T population, ROBUST in all
  24 runs.

## Known limitations of this benchmark

- Case 6 has no ground truth; its expectation is a judgment (same PBMCs, 6 h treatment), not a
  planted fact, and the paired test now calls a small consistent difference (note B).
- One run of each benchmark; runs are deterministic (seeded), so repetition would reproduce
  these numbers rather than estimate variance.
- Nothing here measures false-negative rates for composition claims on real data, because Kang
  has no known compositional effect to recover.
- The planted time effect is one construction with 4 vs 4 donors, one unit above the floor at
  which a permutation test can reach 0.05; its amplitudes are documented choices, not estimates
  of how large real collection-time effects are, and the clock-gene phases disagree with the
  GTEx day/night table for 3 of 14 genes.
- The Mode B rule (the clock-state index decides) was settled after the true IFN-β claims had
  been seen to move clock genes uniformly; the zero false-alarm count for Mode B should be read
  with that in mind.
- The planted doublet block defeats both Scrublet and the mixture diagnostic; the case is a
  miss, and the tool's correct NOT SUPPORTED on it comes from the population's fragility across
  normalisations, not from doublet detection.

## The first run (v0.1), in brief

Sixteen cases, 10 of 16 by the pre-registered rules. The five true Kang claims were FRAGILE
(15 to 16 of 21 runs), all because `min_genes=500` removes 44 % of the cells of a shallow
dataset and, for the monocytes, pushed the stimulated cluster's whole-dataset marker AUROC from
0.83 to 0.79; the null composition case was a false alarm because the matching found the
monocytes in one condition only in those same runs (0.1 % vs 35 %); the planted doublet
population was ROBUST in all 23 runs because an additive four-marker score cannot distinguish a
monocyte+B doublet from "monocyte or B cell", and Scrublet flagged none of the injected cells;
the four false claims, the over-clustering, QC-deletion and batch cases and both controls came
out as expected. Scrublet against demuxlet without the fallback: AUROC 0.87, recall 0.27,
nothing called in the stimulated lane. The full table is in `results_v0.1_table.md`.
