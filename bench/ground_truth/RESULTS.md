# Ground truth from matched bulk RNA-seq: results

Recluster was scored against bulk RNA-seq of the same animals. Each claim is "gene X changes in
bone-marrow phagocytes 4 h after LPS". The reference is a bulk DESeq2 test on three animals per
condition. It is called the truth below because the plan defines it so, but it is itself an
estimate and may not cover the same cells as the single-cell pool (sections 4 and 6). Everything
was set out in [PLAN.md](PLAN.md) before any single-cell test was run on these data. The 40
claims ran on 2026-11-13 (04:38 to 08:24 UTC) with Recluster 0.6.0 (the version recorded in each
`provenance.json`), the default grid and thresholds, and 14 worker processes on a shared laptop.

Every number not marked † comes from `out/scores.md`, `out/scores.json`, `out/summary.json`,
`claims/claims.json`, `bulk_truth/*.csv` or the per-claim `out/<folder>/provenance.json`.
Numbers marked † come from a post hoc analysis run after scoring. They are computed by
`posthoc.py`, which writes them to `posthoc/numbers.json` with supporting tables in `posthoc/`.
That analysis reproduces the bulk truth and Recluster's pseudobulk values exactly.

## 1. Summary

| caller | sensitivity, pooled | specificity, pooled | specificity, mouse | specificity, rat |
|---|---|---|---|---|
| Recluster (ROBUST) | 20/20 (100 %, 84–100) | 17/20 (85 %, 64–95) | 10/10 (100 %, 72–100) | 7/10 (70 %, 40–89) |
| one pseudobulk analysis (baseline run) | 20/20 (100 %, 84–100) | 17/20 (85 %, 64–95) | 10/10 (100 %, 72–100) | 7/10 (70 %, 40–89) |
| naive cell-level test (baseline run) | 20/20 (100 %, 84–100) | 0/20 (0 %, 0–16) | 0/10 (0 %, 0–28) | 0/10 (0 %, 0–28) |

Brackets give the rate and its Wilson 95 % interval. Sensitivity is 10/10 (100 %, 72–100) in
each species for every caller. No claim was FRAGILE, so counting ROBUST or FRAGILE as a call gives
exactly the ROBUST rates.

- **Recluster called all 20 bulk changes.** This says little: they were drawn from genes with
  |bulk log2FC| ≥ 1 at padj < 0.01, and the naive cell-level test also called all 20.
- **Recluster made 3 false calls among the 20 bulk non-changes**, all in rat: Alkbh7, Ezh2 and
  Son. Each was ROBUST in 21 of 21 runs, so the grid gave no warning on any of them.
- **The naive cell-level test (Wilcoxon on cells) called all 20 non-changes.**
- **Recluster's verdicts equalled its baseline run alone on all 40 claims.** The "one pseudobulk
  analysis" is that baseline run: one paired pseudobulk test on the same cells, after the same
  QC and doublet removal, so it is not an independent analysis. On these genes the grid added no
  accuracy over it. Every claim held in all runs or in none, because in these data almost none of
  the grid's choices change the cells tested (section 5).
- **Post hoc†, the false-call rate is not small.** Applied to every eligible gene, the same paired
  one-gene pseudobulk test called 451 of 1,830 rat (24.6 %) and 158 of 1,752 mouse (9.0 %) bulk
  no-change genes at p < 0.05. Recluster chose the paired design itself; rerun unpaired, none of
  the three false calls is made (section 4).

## 2. What was tested

The design is in [PLAN.md](PLAN.md). This section gives the essentials and the details the plan
leaves out.

The data are from Hagai et al. 2018 (Nature 563:197): bone-marrow mononuclear phagocytes from
three mice and three rats, stimulated with LPS for 4 hours or left unstimulated. Bulk and
single-cell libraries come from the same animals (`mouse1-3`, `rat1-3` in the ArrayExpress SDRFs
of E-MTAB-6754 and E-MTAB-6773); of the datasets Squair et al. 2021 use as ground truth, these are
the only ones where that holds. The single-cell data are the authors' per-animal UMI tables
(`<species><n>_<condition>_filtered_by_cell_cluster0.txt.gz`, E-MTAB-6754), which hold the cells
the authors kept as phagocytes: 17,776 mouse cells (9,195 stimulated, 8,581 unstimulated) and
13,277 rat cells (5,952 and 7,325). The bulk data are the count matrices deposited by Squair et
al. (Zenodo 10.5281/zenodo.7659806, `bulk_rnaseq.tar.gz`).

The bulk truth is pydeseq2 0.5.4 with design `~condition`, stimulated against unstimulated,
default settings. It is unpaired because the deposited matrices number the replicates 1–3 within
each condition without saying which animal each column is. A true change has bulk padj < 0.01
and |log2FC| ≥ 1, and is claimed in the direction of the bulk sign. A no-change gene has bulk
padj ≥ 0.5 and |log2FC| < 0.5. Genes in between are not drawn.

A gene is eligible, by rules that use no differential test, when it is in both the single-cell
table and the bulk matrix under the same symbol, has at least 10 bulk counts in total, and is
detected in at least 10 % of the comparison's cells, both conditions pooled. Ensembl IDs are
mapped with BioMart release 116. A single-cell ID that maps to no symbol, or to a symbol shared
with another ID, is excluded, and so is any symbol that occurs more than once in the bulk matrix
(`prepare.py` drops every copy). The eligible pools were 1,175 change and 1,752 no-change genes
in mouse LPS, and 654 and 1,830 in rat LPS. Each comparison has its own numpy `default_rng(0)`,
which draws 10 change genes and then 10 no-change genes, each from the alphabetically sorted pool.

Each claim is `differential_expression` on population `label_value: phagocyte` (every cell),
condition `stim` against `unst`, `sample: [animal, condition]`, with the default grid and
thresholds. True-change genes are claimed in the bulk direction. No-change genes are claimed with
direction `either`: the claim is that the gene changes, which the bulk says is false. Recluster
chose the paired design `~pair + condition` automatically, because the same three animals appear
in both conditions.

Recluster calls a claim when its verdict is ROBUST. The naive cell-level test is the baseline
run's Wilcoxon test on cells, called when its adjusted p (`cell_padj`) is below 0.05 in the
claimed direction. One pseudobulk analysis is the baseline run's replicate-level test alone,
called the same way; for a one-gene claim its adjusted p equals the raw p. Sensitivity is the
share of true-change claims called and specificity the share of no-change claims not called.

The plan has three dated changes, all recorded on 2026-11-12 before any single-cell test:

1. Some pools had fewer eligible genes than 10, and `prepare.py` draws all of them. Mouse
   poly(I:C) had 7 eligible true-change genes and rat poly(I:C) had 9, which gives 76 claims
   instead of 80.
2. The poly(I:C) bulk truth is weak. Unpaired bulk DESeq2 finds 43 (mouse) and 36 (rat) genes at
   padj < 0.05 for poly(I:C), against 4,965 and 4,145 for LPS, so some poly(I:C) "no change"
   genes may change.
3. Only the LPS comparisons were run, for two reasons: on the shared machine one claim had taken
   about 16 minutes, not the 4 estimated, and the poly(I:C) bulk truth is weak (change 2). That
   leaves 40 claims. The 36 poly(I:C) claim files remain in `claims/` but were neither run nor
   scored.

The time reason did not hold. The 40 LPS claims took 190.8 to 536.1 s each (mean 339.0 s,
`out/summary.json`), not about 16 minutes, so at that rate the 36 poly(I:C) claims would have
added roughly 3.4 hours. The remaining reason is the weak bulk truth. The poly(I:C) claims have
not been run; running them now, after the LPS results are known, would be outside the plan.

## 3. Results

| comparison | caller | sensitivity (95 % CI) | specificity (95 % CI) |
|---|---|---|---|
| mouse LPS | Recluster ROBUST | 10/10 (100 %, 72–100) | 10/10 (100 %, 72–100) |
| mouse LPS | one pseudobulk analysis | 10/10 (100 %, 72–100) | 10/10 (100 %, 72–100) |
| mouse LPS | naive cell-level test | 10/10 (100 %, 72–100) | 0/10 (0 %, 0–28) |
| rat LPS | Recluster ROBUST | 10/10 (100 %, 72–100) | 7/10 (70 %, 40–89) |
| rat LPS | one pseudobulk analysis | 10/10 (100 %, 72–100) | 7/10 (70 %, 40–89) |
| rat LPS | naive cell-level test | 10/10 (100 %, 72–100) | 0/10 (0 %, 0–28) |

Pooled rates are in section 1.

**Misses and false calls**

- **Recluster and one pseudobulk analysis:** no misses. False calls: rat Alkbh7, Ezh2 and Son.
- **Naive cell-level test:** no misses. False calls: all 20 no-change genes.

### Every claim

"Pseudobulk" columns are the baseline run's replicate-level test. "Cell" columns are its
Wilcoxon test on cells. A cell padj of 0 is below the floating-point range. "Detected" is the
share of cells with at least one UMI in the baseline run.

| claim | bulk class | bulk log2FC | bulk padj | verdict | runs held | pseudobulk log2FC | pseudobulk p | cell log2FC | cell padj | detected stim / unst | Recluster | one pseudobulk | naive cell |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mouse_lps4/Acsl1 | change | 4.23 | 3e-142 | ROBUST | 23/23 | 4.48 | 2.2e-151 | 1.972 | 0 | 88.1 / 17.4 % | call | call | call |
| mouse_lps4/Ampd3 | change | 3.01 | 1.4e-21 | ROBUST | 23/23 | 3.16 | 4.1e-34 | 0.574 | 0 | 42.3 / 7.1 % | call | call | call |
| mouse_lps4/Azi2 | change | 1.74 | 7.2e-09 | ROBUST | 23/23 | 1.85 | 8.1e-20 | 0.611 | 0 | 70.7 / 34.7 % | call | call | call |
| mouse_lps4/Cry1 | change | 1.25 | 5.1e-05 | ROBUST | 23/23 | 1.27 | 8.2e-08 | 0.137 | 0 | 22.7 / 11.2 % | call | call | call |
| mouse_lps4/Flvcr1 | change | -1.95 | 1.6e-06 | ROBUST | 23/23 | -2.05 | 9.2e-23 | -0.794 | 0 | 14.6 / 49.6 % | call | call | call |
| mouse_lps4/Gprc5b | change | -2.23 | 7.6e-08 | ROBUST | 23/23 | -2.30 | 1.2e-19 | -0.437 | 0 | 5.0 / 22.6 % | call | call | call |
| mouse_lps4/Mtmr14 | change | 2.31 | 1.9e-20 | ROBUST | 23/23 | 2.31 | 2e-25 | 0.337 | 0 | 35.1 / 10.6 % | call | call | call |
| mouse_lps4/Ppfia3 | change | 4.69 | 4.3e-44 | ROBUST | 23/23 | 5.37 | 1.1e-75 | 0.964 | 0 | 54.5 / 2.5 % | call | call | call |
| mouse_lps4/Spred1 | change | 3.40 | 2.4e-37 | ROBUST | 23/23 | 3.53 | 9.8e-54 | 1.124 | 0 | 73.7 / 13.6 % | call | call | call |
| mouse_lps4/Tap1 | change | 3.36 | 1.5e-15 | ROBUST | 23/23 | 3.55 | 3.3e-82 | 0.843 | 0 | 57.4 / 8.9 % | call | call | call |
| mouse_lps4/2210016L21Rik | no change | -0.09 | 0.84 | NOT SUPPORTED | 0/23 | -0.14 | 0.52 | -0.079 | 1.8e-07 | 15.1 / 17.3 % | - | - | **false call** |
| mouse_lps4/Anapc5 | no change | -0.16 | 0.66 | NOT SUPPORTED | 0/23 | -0.21 | 0.29 | -0.270 | 0 | 36.6 / 44.1 % | - | - | **false call** |
| mouse_lps4/Fam193a | no change | -0.16 | 0.71 | NOT SUPPORTED | 0/23 | -0.03 | 0.91 | -0.061 | 5e-05 | 11.1 / 12.7 % | - | - | **false call** |
| mouse_lps4/Itgb7 | no change | -0.28 | 0.56 | NOT SUPPORTED | 0/23 | -0.24 | 0.22 | -0.161 | 9.8e-12 | 18.3 / 21.3 % | - | - | **false call** |
| mouse_lps4/Nmt1 | no change | 0.06 | 0.86 | NOT SUPPORTED | 0/23 | 0.19 | 0.27 | -0.149 | 4e-08 | 51.8 / 50.7 % | - | - | **false call** |
| mouse_lps4/Rab14 | no change | 0.10 | 0.69 | NOT SUPPORTED | 0/23 | 0.22 | 0.19 | -0.187 | 2.8e-38 | 81.8 / 79.6 % | - | - | **false call** |
| mouse_lps4/Selenok | no change | 0.16 | 0.63 | NOT SUPPORTED | 0/23 | 0.25 | 0.17 | -0.203 | 0 | 91.0 / 88.5 % | - | - | **false call** |
| mouse_lps4/Snrpa | no change | -0.21 | 0.66 | NOT SUPPORTED | 0/23 | -0.12 | 0.58 | -0.097 | 1.2e-07 | 17.2 / 19.3 % | - | - | **false call** |
| mouse_lps4/Tbl1x | no change | 0.05 | 0.87 | NOT SUPPORTED | 0/23 | 0.06 | 0.80 | -0.094 | 3e-06 | 17.5 / 19.3 % | - | - | **false call** |
| mouse_lps4/Uimc1 | no change | 0.01 | 0.97 | NOT SUPPORTED | 0/23 | 0.06 | 0.75 | -0.093 | 0.00096 | 21.0 / 21.6 % | - | - | **false call** |
| rat_lps4/Ada | change | 1.53 | 1.2e-11 | ROBUST | 21/21 | 1.72 | 1.3e-36 | 0.324 | 0 | 65.1 / 42.1 % | call | call | call |
| rat_lps4/Ankrd33b | change | 1.43 | 2.9e-10 | ROBUST | 21/21 | 1.60 | 1.8e-25 | 0.203 | 0 | 46.8 / 28.9 % | call | call | call |
| rat_lps4/Baz1b | change | -1.01 | 3.6e-09 | ROBUST | 21/21 | -1.09 | 3.3e-13 | -0.227 | 0 | 21.1 / 43.2 % | call | call | call |
| rat_lps4/Cmpk2 | change | 6.60 | 1.1e-09 | ROBUST | 21/21 | 8.46 | 1.3e-38 | 2.314 | 0 | 71.9 / 2.8 % | call | call | call |
| rat_lps4/Etv1 | change | -1.89 | 2.2e-16 | ROBUST | 21/21 | -2.10 | 3e-38 | -0.375 | 0 | 13.2 / 47.1 % | call | call | call |
| rat_lps4/Fnip2 | change | -1.76 | 1e-22 | ROBUST | 21/21 | -1.96 | 1.4e-19 | -0.131 | 0 | 5.0 / 19.8 % | call | call | call |
| rat_lps4/Mtmr10 | change | -1.03 | 2.9e-05 | ROBUST | 21/21 | -1.44 | 2.7e-16 | -0.133 | 0 | 8.4 / 24.1 % | call | call | call |
| rat_lps4/Polr1e | change | 2.53 | 1.2e-35 | ROBUST | 21/21 | 2.57 | 1.2e-54 | 0.153 | 0 | 29.9 / 8.0 % | call | call | call |
| rat_lps4/Sp110 | change | 1.06 | 0.0036 | ROBUST | 21/21 | 1.79 | 1.3e-29 | 0.606 | 0 | 87.2 / 62.9 % | call | call | call |
| rat_lps4/Tapbp | change | 1.37 | 0.0012 | ROBUST | 21/21 | 1.74 | 1.4e-28 | 0.180 | 0 | 41.7 / 20.5 % | call | call | call |
| rat_lps4/Aarsd1 | no change | -0.08 | 0.70 | NOT SUPPORTED | 0/21 | -0.19 | 0.16 | -0.193 | 0 | 45.7 / 57.7 % | - | - | **false call** |
| rat_lps4/Alkbh7 | no change | 0.10 | 0.91 | ROBUST | 21/21 | 0.35 | 0.012 | -0.039 | 2.1e-06 | 36.7 / 38.5 % | **false call** | **false call** | **false call** |
| rat_lps4/Ezh2 | no change | -0.32 | 0.57 | ROBUST | 21/21 | -0.68 | 2.2e-05 | -0.127 | 3.8e-37 | 16.3 / 24.1 % | **false call** | **false call** | **false call** |
| rat_lps4/Ing3 | no change | 0.13 | 0.69 | NOT SUPPORTED | 0/21 | 0.10 | 0.60 | -0.027 | 4.9e-05 | 13.2 / 15.2 % | - | - | **false call** |
| rat_lps4/Nipsnap3b | no change | -0.04 | 0.88 | NOT SUPPORTED | 0/21 | 0.02 | 0.87 | -0.230 | 0 | 74.1 / 80.6 % | - | - | **false call** |
| rat_lps4/Prpf6 | no change | 0.08 | 0.62 | NOT SUPPORTED | 0/21 | 0.18 | 0.19 | -0.071 | 5.7e-14 | 45.8 / 48.4 % | - | - | **false call** |
| rat_lps4/Serinc2 | no change | -0.04 | 0.87 | NOT SUPPORTED | 0/21 | -0.03 | 0.84 | -0.103 | 6e-22 | 32.5 / 38.0 % | - | - | **false call** |
| rat_lps4/Son | no change | -0.03 | 0.84 | ROBUST | 21/21 | 0.32 | 0.014 | -0.149 | 0 | 93.2 / 94.8 % | **false call** | **false call** | **false call** |
| rat_lps4/Tbl3 | no change | 0.09 | 0.65 | NOT SUPPORTED | 0/21 | 0.17 | 0.31 | -0.032 | 0.00033 | 17.6 / 19.2 % | - | - | **false call** |
| rat_lps4/Ubap2 | no change | 0.05 | 0.81 | NOT SUPPORTED | 0/21 | -0.03 | 0.84 | -0.079 | 1.9e-21 | 26.9 / 32.5 % | - | - | **false call** |

**The naive cell-level test** called every no-change gene, and every one has a negative cell
log2FC: from −0.270 to −0.061 in mouse and from −0.230 to −0.027 in rat. The largest cell padj
among them is 0.00096 (mouse Uimc1); in rat it is 0.00033 (Tbl3). Detection is lower in
stimulated cells for 17 of the 20 genes, including all 10 in rat.

The cell-level test treats each cell as a replicate, so with about 6,000 to 9,000 cells per
group, shifts of 0.027 to 0.270 log2 reach significance. Why all 20 shifts are negative is not
tested here. One possibility is the library-size normalisation of the cell-level test: in rat,
stimulated cells have larger libraries (median 24,452 against 21,477 UMI†). The smallest |cell
log2FC| among the true changes is 0.131 (rat Fnip2), so effect size does not separate the two
classes in this test either.

## 4. The three false calls

Bulk values come from `bulk_truth/rat_lps4.csv` and single-cell values from
`out/rat-lps4-*/provenance.json`. Entries marked † come from the post hoc analysis. That analysis
reproduces the bulk DESeq2 exactly, and Recluster's paired pseudobulk on all 13,277 rat cells
exactly (for example, Alkbh7 log2FC 0.3545, p 0.01189, genome-wide padj 0.0293).

| | Alkbh7 | Ezh2 | Son |
|---|---|---|---|
| bulk log2FC (lfcSE†) | 0.0958 (0.471) | -0.3152 (0.372) | -0.0299 (0.088) |
| bulk baseMean / padj / total counts | 259.18 / 0.9103 / 1,574 | 192.94 / 0.5742 / 1,223 | 8,439.76 / 0.8435 / 53,922 |
| pseudobulk log2FC, baseline run | 0.354338 | -0.681578 | 0.317053 |
| pseudobulk p (claim padj), baseline | 0.0118716 | 2.21e-05 | 0.013611 |
| pseudobulk padj, genome-wide, baseline | 0.029218 | 9.8e-05 | 0.032805 |
| range across the 21 runs: log2FC | 0.3543 to 0.3553 | -0.6910 to -0.6816 | 0.3171 to 0.3179 |
| range across the 21 runs: claim padj | 0.01169 to 0.01189 | 1.79e-05 to 2.25e-05 | 0.01346 to 0.01361 |
| detected, stim / unst | 36.7 % / 38.5 % | 16.3 % / 24.1 % | 93.2 % / 94.8 % |
| pseudobulk sign against bulk sign | same (+) | same (−) | opposite (bulk −0.03, single cell +0.32) |
| per-animal pseudobulk log2(stim/unst), DESeq2 size factors† (rat1 / rat2 / rat3) | 0.173 / 0.479 / 0.418 | -0.506 / -0.725 / -0.860 | 0.177 / 0.541 / 0.235 |
| per-animal log2(stim/unst), library-size (CP10k) scaling† | -0.373 / -0.073 / -0.161 | -1.051 / -1.277 / -1.438 | -0.369 / -0.011 / -0.343 |

No grid choice moves any of the three calls: across the 21 runs the claim padj stays between
0.0117 and 0.0136 for Alkbh7 and Son, and at or below 2.25e-05 for Ezh2. With the baseline run's
genome-wide correction each would also pass at padj 0.05.

### What the data show

**The calls depend on the paired design.**

- Recluster's pseudobulk test pairs the stimulated and unstimulated samples of each animal
  (`~pair + condition`). The bulk truth is unpaired.
- Rerun unpaired on the all-cells pseudobulk counts (before doublet removal)†, none of the three
  would be called: p 0.326 (Alkbh7), 0.229 (Ezh2) and 0.307 (Son).
- No grid axis touches the pseudobulk design or its normalisation, so the grid could not show
  this.

**For Son, single cell and bulk disagree.**

- The bulk estimate is precise: lfcSE 0.088, 95 % interval about −0.20 to 0.14. That interval
  excludes the single-cell +0.32.
- Pairing the bulk columns by replicate index† does not change this: −0.029 (SE 0.100), p 0.77.
  This assumes that index i is the same animal in both conditions. Which animal each bulk column
  is is not recorded, and the plan declined to assume it.
- Son is not linked to cell cycle (cycling association −0.10 log2†).

**For Alkbh7, bulk is too imprecise to rule out the single-cell estimate.**

- Unpaired bulk lfcSE is 0.471. Bulk paired by index† gives 0.138 (SE 0.193), p 0.47. Both
  intervals contain +0.35.
- Alkbh7 is not linked to cell cycle (association 0.03†).

**For Ezh2, the fall is consistent with a loss of cycling cells.**

- On a post hoc definition, a cell counts as cycling when at least 4 of 17 G2M genes or 4 of 14
  S-phase genes are detected† (the list has 18 G2M genes, but Ccnb2 is not in the rat single-cell
  table). By that definition the share of cycling rat cells falls after LPS†:

  | | rat1 | rat2 | rat3 |
  |---|---|---|---|
  | unstimulated | 56.9 % | 29.0 % | 34.3 % |
  | stimulated | 16.6 % | 9.5 % | 11.1 % |

  The definition depends on sequencing depth. Mouse libraries are shallower (median 9,761
  against 22,611 UMI per cell in rat†), and on the same rule, with all 18 G2M genes, 0.5–1.7 % of
  mouse cells count as cycling†. That figure is not comparable with rat.
- In unstimulated rat cells, Ezh2 is 5.6 times (log2 2.48) higher in cycling cells† (6.4 times
  without the pseudocount†). A gene's cycling association, here and above, is
  log2((CP10k in cycling + 0.01) / (CP10k in non-cycling + 0.01)) over unstimulated cells; the
  pseudocount of 0.01 CP10k keeps genes absent from one group finite. Within
  non-cycling cells, with library-size (CP10k) scaling, its per-animal log2(stim/unst) is +0.672,
  +0.004 and +0.046†: two animals show no fall and rat1 shows a rise. On the same scaling, all
  cells give −1.051, −1.277 and −1.438 (table above).
- Bulk shows a much smaller fall in cell-cycle genes than the pseudobulk†: Top2a −0.90 (bulk padj
  0.040) against −4.43, Mki67 −0.99 (padj 0.089) against −3.71, and Cdk1, Ube2c and Ccnb1 −0.36,
  −0.20 and −0.21 against −2.90, −2.85 and −3.27. The last three are in the bulk no-change class.
- The bulk Ezh2 estimate (−0.32, lfcSE 0.372) has the same sign as the single cell but is far
  from significant (padj 0.57). Bulk paired by replicate index, an animal order the plan declined
  to assume, gives −0.360 (SE 0.214), p 0.092†.
- Under the plan's rule Ezh2 is a false call.

**The sign for Alkbh7 and Son depends on normalisation.**

- This is visible in the committed files: the cell-level test, on library-size-normalised
  (CP10k) data, gives −0.039 (Alkbh7) and −0.149 (Son); the pseudobulk test, with DESeq2 size
  factors, gives +0.354 and +0.317.
- Per animal, all three animals are positive for both genes with DESeq2 size factors and negative
  with library-size scaling†.
- The two scalings differ by a per-animal offset of 0.546, 0.552 and 0.578 log2 in rat†, and of
  0.40–0.67 in mouse†. Stimulated rat cells have larger libraries (median 24,452 against 21,477
  UMI) but fewer genes detected (3,211 against 3,448)†.
- The bulk truth also uses DESeq2 size factors, so this does not explain why bulk and pseudobulk
  disagree for Son. It shows only that the single-cell estimate is close to zero and that its
  sign rests on how stimulated and unstimulated libraries are scaled.

### What these files cannot decide

- **Whether Alkbh7 and Son really change in the cells.** Bulk and pseudobulk use the same
  normalisation method yet disagree for Son; for Alkbh7 the bulk is too imprecise to say.
- **Why the bulk shows a much smaller fall in cell-cycle genes.** The single-cell files are the
  authors' `filtered_by_cell_cluster0` tables, so one cluster was kept. How many cells were
  excluded cannot be recovered. The SDRF of E-MTAB-6773 lists two bulk libraries per animal and
  condition (for example rat1_LPS4 and rat1_LPS4a). The Squair matrix has three columns per
  condition and does not say how the libraries were combined. The bulk sample may therefore not
  be the same cell population as the single-cell pool. That is consistent with the cycling
  discrepancy but cannot be measured here.
- **Ambient RNA.** Only filtered cells are available, with no empty droplets or raw matrices, so
  it cannot be tested.

### How typical the three are (post hoc, †)

Recluster's one-gene paired pseudobulk test (p < 0.05) was applied, on all cells before doublet
removal, to every eligible gene, not only the 40 drawn:

| | no-change genes called, paired | no-change genes called, unpaired | change genes called with bulk sign, paired | change genes called with bulk sign, unpaired |
|---|---|---|---|---|
| rat | 451/1,830 (24.6 %) | 70/1,830 (3.8 %) | 653/654 | 646/654 |
| mouse | 158/1,752 (9.0 %) | 69/1,752 (3.9 %) | 1,171/1,175 | 1,150/1,175 |

- **The draw matched the population rates.** At these rates, a draw of 10 is expected to give
  about 2.5 false calls in rat and 0.9 in mouse. The observed numbers were 3 and 0.
- **Bulk noise alone does not explain the rat rate.** For rat no-change genes, the median lfcSE
  is 0.167 for unpaired bulk and 0.142 for paired pseudobulk. If the paired test were calibrated
  and these genes did not change in the single-cell pool, it would call about 5 %, not 24.6 %.
  These files do not separate two explanations: many genes differ between the single-cell pool
  and the bulk by 0.2–0.5 log2, or the paired test (3 pairs, so 2 residual degrees of freedom) is
  anti-conservative here. Mouse, at 9.0 %, is also above 5 %.
- **The cycling-cell loss explains part of it.** Across the 1,830 rat no-change genes,
  enrichment in cycling cells correlates with the gap between pseudobulk and bulk log2FC
  (Spearman −0.404, p 1.2e-72). Of the 58 genes at least 2× higher in cycling cells, 48 are
  called, all of them down. In mouse the same correlation is −0.046.
- **It does not explain all of it.** Among rat genes not linked to cycling (|association| <
  0.25), 288 of 1,296 are still called (194 up, 94 down). Their median pseudobulk-minus-bulk
  shift is +0.207; across all 1,296 it is +0.043.
- **The unpaired test is not simply better.** The unpaired pseudobulk calls fewer no-change
  genes, but its median lfcSE is 0.312 in rat, so it is less powerful on everything. These
  numbers are post hoc and are not a rule to adopt.
- **Effect size separates the classes only partly.** Alkbh7 and Son have single-cell |log2FC|
  below 0.5, inside the bulk "no change" band. Ezh2 (0.68) does not. All 20 true-change claims
  have pseudobulk |log2FC| ≥ 1.09, but the truth rule |bulk log2FC| ≥ 1 nearly guarantees that.

## 5. What the grid added

On all 40 claims the grid verdict equals the baseline pseudobulk call. Every claim held in 23/23
or 0/23 runs (mouse), or 21/21 or 0/21 (rat), and the largest spread of pseudobulk log2FC across
runs on any claim is 0.068 (mouse Gprc5b, −2.299 to −2.231).

The reason is that almost no grid choice changes the cells tested. The population is one label
covering every cell, so clustering, highly variable genes, PCs, resolution, normalisation for
clustering and seed never change them. The authors' tables are already filtered (the fewest genes
in any cell is 1,100 in mouse and 2,425 in rat†), so QC removes no cells except in one setting:
`n_removed_qc` is 0 in every run except mouse's `max_mito_pct=5` runs, which remove 1,130 cells.
Doublet removal is the only other change to the cell set; in the baseline it removes 436 mouse
cells and 210 rat cells. The data have no batch column, so batch correction was off in every run.

Rat also ran without mitochondrial QC. Rat mitochondrial genes are named `Mt-`, and
`detect_mito_prefix` (`recluster/data.py`) tries only `MT-` and `mt-`. Every rat run warns "no
mitochondrial genes found; max_mito_pct has no effect", and the rat grid has 21 runs against
mouse's 23. The highest `Mt-` share in any rat cell is 9.99 %†, so the baseline 10 % cut would
have removed nothing anyway.

The choices that mattered for the false calls are paired against unpaired design, pseudobulk
normalisation and the definition of the population. None of them is a grid axis.

## 6. Limitations

- **Only 40 claims.** The Wilson intervals are wide: specificity 64–95 % pooled and 40–89 % in
  rat.
- **The genes come from the extremes.** They are strong bulk changes (|log2FC| ≥ 1, padj < 0.01)
  against near-zero non-changes (|log2FC| < 0.5, padj ≥ 0.5). The study says little about
  modest effects, where calls are harder.
- **One cell population.** Every claim tested all cells, so population matching, the part of
  Recluster that varies most across the grid, never mattered here.
- **The bulk truth is unpaired.** Recluster's test is paired. The bulk replicates cannot be
  assigned to animals with certainty.
- **Rat ran without mitochondrial QC**, because of the `Mt-` prefix (section 5).
- **The time-of-day check did not run.** It was NOT ASSESSABLE on all 40 claims, with the
  reason "0 of 13 panel genes have at least 10 counts". The real cause is gene names. The clock
  panel is matched by upper-case human symbols (`recluster/time_confound.py`, `panel_expression`),
  while mouse and rat symbols are title case. Twelve of the 13 panel genes are present in both
  files (`Per1`, `Cry1`, `Dbp` and others). If the names had matched, the check could have
  reported a possible time confound and capped a ROBUST verdict at FRAGILE. These files cannot
  say whether it would have. As written, the check never runs on mouse or rat gene names.
- **poly(I:C) was not scored.** Its bulk truth is weak (43 and 36 genes at padj < 0.05). Time was
  also given as a reason, but the LPS claims ran faster than feared (section 2), so time alone
  would not have ruled it out.
- **The bulk truth is itself an estimate.** It comes from three animals per condition, it may
  not cover the same cell population as the single-cell pool, and for genes such as Alkbh7 its
  "no change" rests on wide intervals.
- **The post hoc analysis was defined after scoring.** `posthoc.py` reruns it from the
  committed files and the input data, but its choices (the cycling rule, the 2× and 0.25 cuts,
  the 0.01 CP10k pseudocount) were made after the results were seen, so the † numbers are
  exploratory. Its check against Recluster's own pseudobulk runs only when the git-ignored
  `out/<folder>/provenance.json` files are present.

## 7. How to reproduce

There is no data README for this study; the inputs are as in PLAN.md and `prepare.py`. All
paths are relative to the repository root, and nothing under `data/` is in git.

- `data/hagai2018/sc/<species><n>_<lps4|pic4|unst>_filtered_by_cell_cluster0.txt.gz`: the 18
  per-animal single-cell tables from ArrayExpress E-MTAB-6754.
- `data/hagai2018/bulk_rnaseq.tar.gz` from Zenodo 10.5281/zenodo.7659806, extracted into
  `data/hagai2018/bulk/`. This gives
  `data/hagai2018/bulk/bulk_rnaseq/matrices/Hagai2018_<mouse|rat>-<lps|pic>_{matrix,features,meta}.tsv.gz`.
- `data/hagai2018/biomart/mmusculus.tsv` and `rnorvegicus.tsv`: Ensembl BioMart release 116
  exports with the columns `Gene stable ID` and `Gene name`.
- Optionally, `E-MTAB-6754.sdrf.txt` and `E-MTAB-6773.sdrf.txt`, used to check that bulk and
  single cell come from the same animals. The scripts do not read them.

Then run:

```bash
uv run python bench/ground_truth/prepare.py      # h5ad files, bulk_truth/*.csv, claims/
uv run python bench/ground_truth/run_claims.py   # 40 LPS claims, 14 workers by default; resumable
uv run python bench/ground_truth/score.py        # out/scores.json, out/scores.md
uv run python bench/ground_truth/posthoc.py      # posthoc/: the † numbers, about 1 min
```

`prepare.py` builds all four comparisons and writes four `.h5ad` files of 278–374 MB each.
`run_claims.py` runs only the LPS claims (change 3).

The 40 claims took 13,561.1 s in total (about 3.8 h): 7,029.1 s for mouse and 6,532.0 s for rat,
with 190.8 to 536.1 s per claim (`out/summary.json`). They ran one at a time, each with 14
worker processes, on a laptop that was running other jobs.
