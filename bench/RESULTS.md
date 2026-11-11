# Benchmark results

Measured on 2026-11-11 with Recluster 0.2.0, default grid and default thresholds, on an Apple
M5 Pro laptop (`recluster bench --jobs 4`, 35 minutes wall time). Every number below is copied
from `results.json`; nothing was tuned after seeing it. Where the outcome is a miss or a false
alarm, it says so.

Verdict rule for all cases: ROBUST if the claim held in >= 80 % of grid runs, FRAGILE if in
>= 20 %, otherwise NOT SUPPORTED. "Held" counts runs where the replicate-aware test passed.
"Cell-level" is the share of runs where the naive cell-level statistic (Wilcoxon on cells, or
chi-square on pooled cells) would have called the claim; it never enters the verdict.

## Headline table

| # | Case | Expected | Verdict | Held | Cell-level | Outcome |
|---|---|---|---|---|---|---|
| 1 | Kang: ISG15 up in CD14+ monocytes, stim vs ctrl | ROBUST | FRAGILE | 16/21 | 100 % | missed (see note A) |
| 2 | Kang: ISG15, IFI6, IFIT3, MX1 all up in CD14+ monocytes | ROBUST | FRAGILE | 16/21 | 100 % | missed (note A) |
| 3 | Kang: ISG15 up in CD4 T cells | ROBUST | FRAGILE | 15/21 | 81 % | missed (note B) |
| 4 | Kang: CD14+ monocytes exist across donors | ROBUST | FRAGILE | 16/21 | 76 % | missed (note A) |
| 5 | Kang: NK cells exist across donors | ROBUST | FRAGILE | 15/21 | 71 % | missed (note A) |
| 6 | Kang: CD14+ monocyte proportion differs stim vs ctrl (no known shift) | NOT SUPPORTED | FRAGILE | 8/21 | 100 % ⚑ | false alarm (note C) |
| 7 | False: ISG15 differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 1/21 | 81 % ⚑ | correct |
| 8 | False: LYZ differs between fake donor groups (ctrl only) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 100 % ⚑ | correct |
| 9 | False: monocyte proportion differs between fake donor groups | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 19 % | correct |
| 10 | False: ISG15 up in fake group A (condition shuffled at donor level) | NOT SUPPORTED | NOT SUPPORTED | 0/21 | 100 % ⚑ | correct |
| 11 | Planted: 150 synthetic doublets as a fake population | not ROBUST, doublet axis breaks it | ROBUST | 23/23 | 100 % | **missed** (note D) |
| 12 | Planted: over-clustering shard of CD4 T cells | not ROBUST | NOT SUPPORTED | 0/23 | 0 % | caught |
| 13 | Control: a shard that is a whole population (B cells) | ROBUST | ROBUST | 23/23 | 100 % | control OK |
| 14 | Planted: QC threshold deletes NK cells | failure table names max_mito_pct | NOT SUPPORTED | 4/23 | 17 % | caught |
| 15 | Planted: batch artifact claimed as a CD4 T subtype | not ROBUST, batch axis breaks it | NOT SUPPORTED | 4/24 | 17 % | caught |
| 16 | Control: the real CD4 T population under the batch effect | ROBUST | ROBUST | 24/24 | 100 % | control OK |

⚑ = pseudoreplication flag raised by the report (cell-level statistics call it in >= 80 % of
runs while replicate-aware statistics call it in < 50 %).

Score by the pre-registered rules: 10 of 16 as expected (4 false claims correct, 3 planted
errors caught, 2 controls OK, 1 more if the null composition case is read as "not ROBUST");
5 true claims came out FRAGILE instead of ROBUST, 1 planted error missed, 1 false alarm.

## Notes on the misses

**A. The five true Kang claims are FRAGILE for one reason: `min_genes=500`.** Every one of the
five failing runs for cases 1, 2, 4 and 5 (and five of the six for case 5) has `min_genes=500`;
all runs with 100 or 200 hold. Kang 2018 is shallow: the median cell has 519 detected genes,
and a 500-gene filter keeps only 56 % of cells, unevenly by type (CD14+ monocytes 77 %, NK 34 %,
CD8 T 23 %, megakaryocytes 10 %). Recluster's Jaccard rule counts cells removed by QC as
absent from the population, so under that filter the matched monocyte set overlaps the
reference by less than 50 % and the claim is scored as not holding, even though the
pseudobulk test on the surviving cells still finds ISG15 up. Read as "the claim's survival
depends on QC" this is a true statement about the dataset; read against the pre-registered
expectation (ROBUST) it is a miss, and it exposes two things worth fixing in v0.3:
the fixed grid value 500 is not a defensible threshold for a dataset whose median is 519
(QC axes should be relative to the data, e.g. quantiles), and existence should probably be
judged on the cells a run kept rather than penalised for cells QC removed. Neither change
was made for this report.

**B. CD4 T cells (case 3)** fail in six runs: `doublet_removal=off`, `normalization=pearson`,
`n_hvgs=1000`, `resolution=0.3`, `resolution=0.5` and one random combination. All are
matching failures, not statistical ones: CD4 T cells are defined by IL7R, CD3D and CCR7, which
separate them from CD8 T cells only weakly, and at low resolution or with fewer features the
cluster that scores highest for those markers also contains CD8 and other T cells, so its
marker AUROC drops below 0.80. Using the labels as given, ISG15 is up in CD4 T cells by both
tests (pseudobulk log2FC 4.9, padj 4e-99).

**C. The null composition claim (case 6) is a false alarm by the pre-registered rule.** With the
labels as given there is no shift (logit t-test p = 0.67; Mann-Whitney p = 0.72; the pooled
cell-level chi-square says p = 0.008, which is the pseudoreplication the benchmark is about).
The 8 runs in which the claim held are the five `min_genes=500` runs plus `normalization=pearson`,
`resolution=2.0` and one random combination. Under `min_genes=500` cell loss is
type-dependent and slightly condition-dependent (52 % of ctrl vs 60 % of stim cells kept), and
the monocyte share among retained cells becomes 37 % in ctrl vs 31 % in stim: a real
compositional difference, created by the QC choice. Recluster reports it as FRAGILE with those
runs listed, which is the right description of what happened, but the claim had no biological
basis, so it is counted as a false alarm.

**D. The planted doublet population was missed outright.** Scrublet flagged 0 of the 150
injected cells at the baseline QC (false-positive rate on real cells 1.0 %), so no grid run
removed them and the fake CD14+ LYZ+ MS4A1+ CD79A+ cluster held in all 23 runs. The
mechanism is known: Scrublet scores a cell by whether its neighbours are simulated doublets,
and 150 near-identical monocyte+B doublets are each other's neighbours. Their median
Scrublet score is about twice the dataset median, which the report now shows as an
informational number, but no verdict rule uses it. A heterotypic-marker test (a population
co-expressing the defining markers of two established populations) is the obvious next check
and is on the roadmap.

## How each case was built

### Real data (Kang et al. 2018, Nat Biotechnol 36:89)

Processed object from pertpy (`https://exampledata.scverse.org/pertpy/kang_2018.h5ad`,
SHA-256 `e6a5adac…9830`): 24,673 cells, 15,706 genes, raw counts; `label` (ctrl 12,315 / stim
12,358), `replicate` (8 donors), `cell_type` (8 types). The object has no mitochondrial genes,
so the `max_mito_pct` axis is dropped and the grid has 21 runs (14 one-at-a-time, 6 random).
Replicate = donor × condition (16 samples). No batch column is declared because condition is
confounded with 10x lane.

- Cases 1 to 5: unmodified data. Populations are located by markers (CD14/LYZ/S100A8;
  IL7R/CD3D/CCR7; GNLY/NKG7/KLRD1 minus CD3D) with `cell_type` as the reference set.
- Case 6: unmodified data; a composition claim with no planted or known effect.
- Cases 7 to 9: control cells only; donors shuffled with seed 0 and split into
  A = patient_1016, patient_1039, patient_107, patient_1256 and
  B = patient_101, patient_1015, patient_1244, patient_1488. Nothing differs between A and B
  by construction.
- Case 10: all cells; each of the 16 donor × condition samples assigned to A or B so that each
  fake group holds 4 ctrl and 4 stim samples (seed 0).

### Doublets vs demuxlet (GEO GSE96583, batch 2)

The processed object holds singlets only, so Scrublet was compared with demuxlet on the two
original lanes (29,065 barcodes; 24,679 singlets, 3,169 doublets, 1,217 ambiguous; 313
barcodes without a call). Scrublet ran exactly as in the pipeline: per lane, baseline QC, exact
neighbours, automatic threshold. Ambiguous barcodes are excluded.

| Lane | Cells | demuxlet doublet rate | Scrublet called | Precision | Recall | AUROC of score |
|---|---|---|---|---|---|---|
| both | 27,406 | 11.5 % | 4.6 % | 0.67 | 0.27 | 0.87 |
| ctrl | 13,853 | 11.5 % | 9.1 % | 0.67 | 0.53 | 0.88 |
| stim | 13,553 | 11.5 % | 0.0 % | – | 0.00 | 0.86 |

The score ranks demuxlet doublets well in both lanes (AUROC 0.86 to 0.88), but Scrublet's
automatic threshold called no doublets at all in the stimulated lane. demuxlet can only see
doublets from two different individuals (about 7/8 of doublets with 8 donors), so precision is
a lower bound. Conclusion as measured: with the automatic threshold, doublet removal in
Recluster is unreliable, and the `doublet_removal` axis should be read accordingly. A threshold
tied to the expected doublet rate is a candidate fix; it was not applied here.

### Planted errors (PBMC3k, `scanpy.datasets.pbmc3k()`, 2,700 cells)

Reference labels come from Recluster's baseline pipeline at resolution 1.0, named by markers
(CD4 T 1,164; CD14 Mono 450; B 332; CD8 T 313; FCGR3A Mono 186; NK 160; DC 31; Mk 15;
49 filtered). Grid: 23 runs (24 when a batch column exists). Seed 0 throughout.

- Case 11: 150 synthetic doublets, each the summed raw counts of a random CD14 monocyte and a
  random B cell, appended as new cells. Claim: a population with markers CD14, LYZ, MS4A1,
  CD79A exists. Caught if not ROBUST and every holding run had doublet removal off.
- Cases 12 and 13: Leiden at resolution 3.0 gives 30 clusters; CD4 T cells are split into 16
  shards. The largest shard (cluster_14, 116 cells, 100 % CD4 T) is claimed as a population via
  labels in audit mode. Control: cluster_2 (328 cells, Jaccard 0.99 with B cells). Caught if
  not ROBUST; control passes if ROBUST.
- Case 14: mitochondrial counts of the 160 reference NK cells rescaled so each cell's mito
  fraction is uniform in 12 to 16 % (all other cells stay below 6 %). The baseline
  `max_mito_pct=10` removes all of them. Claim: NK cells (GNLY, NKG7, KLRD1; not CD3D) exist.
  Caught if the top row of the failure table is `max_mito_pct`. Result: the claim held only in
  the 4 runs with `max_mito_pct=20`; the one-at-a-time flip is `max_mito_pct=20`.
- Cases 15 and 16: a random 40 % of cells assigned to "batch2"; 300 expressed genes multiplied
  by 2.5 and 300 others divided by 2.5 in those cells. Clustering without correction yields
  cluster 1 (450 cells, 98 % batch2, CD4 T), labelled "CD4 T subtype" and claimed in audit mode
  with `batch` declared. Caught if not ROBUST and every holding run had batch correction off
  (the 4 holding runs all had it off; the one-at-a-time flip is `batch_correction=off`).
  Control: the reference CD4 T population, ROBUST in all 24 runs.

## Known limitations

- The grid's absolute QC values are not adapted to the dataset; on a shallow dataset
  `min_genes=500` is an outlier choice, not a defensible alternative, and it drives notes A
  and C.
- Existence is scored with a Jaccard that penalises QC removal, so "QC deletes part of the
  population" and "the population is not there" are not distinguished in the verdict, only in
  the failure table.
- Doublet removal relies on Scrublet's automatic threshold, which failed on one of two real
  lanes and on a block of identical synthetic doublets.
- Marker-based matching of CD4 T cells is weak at low resolution; label-mode matching would
  have been more stable but was not used for the true claims because the point was to find
  the population independently of the labels.
- Case 6 has no ground truth; its expectation is a judgment (same PBMCs, 6 h treatment), not a
  planted fact.
- One run of each benchmark; Kang runs are deterministic (seeded), so repetition would
  reproduce these numbers rather than estimate variance.
- Nothing here measures false-negative rates for composition claims on real data, because
  Kang has no known compositional effect to recover.

## Threshold and design changes made during the benchmark

Recorded for honesty; none changes a number above.

- Before any benchmark claim ran, a probe of the doublet construction showed Scrublet would
  miss the injected block. The report gained an informational median-Scrublet-score line.
  No threshold or rule was changed.
- Scrublet was switched from scanpy's `batch_key` loop to Recluster's own per-sample loop
  with exact neighbours, because the scanpy path crashes with scipy 1.17 on large samples.
- The `max_mito_pct` axis is skipped when the input has no mitochondrial genes.
