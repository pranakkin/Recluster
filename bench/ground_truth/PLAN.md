# Ground truth from matched bulk RNA-seq: plan

Written before any bulk or single-cell differential expression was computed on these data. Once
the claims are drawn, nothing below changes; any later change is added under "Changes" with its
reason and the date.

## Data

Hagai et al. 2018 (Nature 563:197), bone-marrow mononuclear phagocytes stimulated for 4 hours
with LPS or poly(I:C), or left unstimulated. Single-cell and bulk RNA-seq were made from the same
animals (`mouse1-3`, `rat1-3` in the ArrayExpress SDRFs of E-MTAB-6754 and E-MTAB-6773). Squair et
al. 2021 (Nat Commun 12:5692) used these as ground truth for single-cell differential expression.

- Single cell: the authors' per-animal UMI tables, `<species><n>_<condition>_filtered_by_cell_cluster0.txt.gz`
  from E-MTAB-6754 (cells the authors kept as phagocytes).
- Bulk: the count matrices Squair et al. deposited (Zenodo 10.5281/zenodo.7659806,
  `bulk_rnaseq.tar.gz`, `Hagai2018_<species>-<lps|pic>_*`); E-MTAB-6773 has no processed counts.
- Four comparisons, each 3 vs 3 animals: mouse LPS, mouse poly(I:C), rat LPS, rat poly(I:C), each
  against unstimulated cells of the same animals.

Of the datasets Squair et al. call ground truth, these are the only ones whose bulk and single-cell
samples come from the same individuals (checked in the SDRFs). Cano-Gamez 2020 and Angelidis 2019
use different donors or mice for bulk and single cell; Reyfman 2019 provides only a published
bulk table. They are not used.

## Bulk truth

DESeq2 (pydeseq2, the version pinned in uv.lock) on the bulk counts of each comparison, design
`~condition` (unpaired), stimulated vs unstimulated, default settings. Unpaired because the
deposited matrices number the replicates 1-3 within each condition without saying which animal
each column is; pairing them would rest on that assumption. Classes:

- **true change**: bulk padj < 0.01 and |log2 fold change| >= 1; the direction is the sign;
- **no change**: bulk padj >= 0.5 and |log2 fold change| < 0.5;
- genes in between belong to neither class and are not drawn.

## Which genes become claims

Eligible genes, decided from data that do not involve any differential test: the gene is in both
the single-cell table and the bulk matrix (same symbol; the single-cell tables use Ensembl IDs,
mapped to symbols with Ensembl BioMart release 116, and an ID that maps to no symbol or to a symbol
shared with another ID is not eligible), has at least 10 bulk counts in total, and
is detected (UMI > 0) in at least 10 % of the single cells of the comparison, both conditions
pooled. The same rule applies to both classes, so "no change" genes are not just unexpressed ones.

From the eligible genes of each comparison, draw at random (numpy `default_rng(0)`, genes sorted
alphabetically before drawing) N true-change genes and N no-change genes. N is fixed from runtime
before drawing (section "N"), not from any result.

Each drawn gene is one claim spec: `differential_expression`, population = all cells (the files
hold one cell type), condition stimulated vs unstimulated, sample = animal x condition, default
grid and thresholds. True-change genes are claimed in the bulk direction (`up_in_a` or
`down_in_a`); no-change genes are claimed with direction `either` (the claim that the gene
changes, which the bulk says is false).

## Scoring

- Recluster calls a claim when its verdict is ROBUST. Also reported: ROBUST or FRAGILE.
- Naive cell-level testing: the Wilcoxon test on cells in the baseline run (Recluster records it
  as `cell_padj`), called at p < 0.05 in the claimed direction.
- A single pseudobulk analysis: the replicate-level test in the baseline run alone, called the
  same way (shows what the grid adds over one well-chosen analysis).
- Sensitivity = calls among true-change claims; specificity = no call among no-change claims;
  each with a Wilson 95 % interval, per comparison and pooled. Every miss and false call is
  listed.

## N

N = 10 per class per comparison: 80 claims (4 comparisons x 2 classes x 10). Fixed from runtime,
not from a result: each comparison has about 12,000 cells, half of Kang 2018, where one claim takes
about 4 minutes on this laptop, so 80 claims take about 3 hours. The machine was busy with another
run, so this is an estimate rather than a timing run.

## Changes

Recorded after the draw (2026-11-12), before any single-cell test was run:

1. **Fewer eligible genes than N.** The plan did not say what happens when a class has fewer than
   N eligible genes. `prepare.py` draws all of them. This happened for true-change genes in mouse
   poly(I:C) (7 eligible) and rat poly(I:C) (9), so there are 76 claims, not 80.
2. **The poly(I:C) bulk truth is weak.** Unpaired bulk DESeq2 finds 43 (mouse) and 36 (rat) genes
   at padj < 0.05 for poly(I:C), against 4,965 and 4,145 for LPS. The genes it does find are
   expected ones (Ifnb1, Mx2, Ifna5, Egr1, Tnf). So the bulk test is underpowered for poly(I:C),
   and some of its "no change" genes may in fact change. A Recluster call on such a gene would
   count as a false call. Results are reported per comparison as planned; the poly(I:C)
   specificity has to be read with this in mind, and the LPS comparisons are the cleaner test.
3. **Only the LPS comparisons are run (2026-11-12).** With the machine shared, one claim took
   about 16 minutes instead of the 4 estimated under "N", so 76 claims would take most of a day.
   The two poly(I:C) comparisons are dropped, for the reason in change 2: their bulk truth is too
   weak to score against. This was decided before any single-cell test was run on these data. The
   claims are 40 (mouse LPS and rat LPS, 10 true-change and 10 no-change genes each, the genes
   already drawn), and the poly(I:C) claim files stay in `claims/` but are not run or scored.
