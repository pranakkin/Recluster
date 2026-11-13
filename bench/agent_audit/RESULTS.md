# Agent audit: results

Claude Opus 5.5 (`claude-opus-5-5`, through the Anthropic API) was given the Kang 2018
interferon-β PBMC counts, with only a stim/ctrl label and donor ids, and asked for a full
analysis. It ran once. The API has no seed, so the run cannot be replayed exactly. We turned
every checkable claim in its final answer into a Recluster claim spec: 54 specs, 47 primary
and 7 secondary. Secondary specs are literal readings of two group statements, tested in the types
for which the agent's own table has no row. The 11 population specs were then run a second time
with the sample column they had lacked (amendment 9, below); that re-run replaces the first as
the primary population result, and both are reported, so 65 specs ran in all. Each spec ran on
the same 21-run grid. ROBUST means the claim held in at least 80 % of runs, FRAGILE in at least
20 %, otherwise NOT SUPPORTED. Runs that could not test a claim count as failures.

- **Populations (11):** 9 ROBUST, 1 FRAGILE (CD4 T memory, 14/21), 1 NOT SUPPORTED
  (Megakaryocyte/platelet, 1/21). The agent listed that label as less certain, but still gave
  "more platelets (0.29% → 0.56%)" as one of its two strongest trends.
- **Expression (25):** 23 ROBUST, 2 FRAGILE (five down-genes in cDC, 5/21; NELL2 down in NK,
  6/21). 9 of the 23 ROBUST specs are the same 17 interferon genes in 9 types. Both FRAGILE
  specs hold when Recluster tests them once on the agent's own labels. With the populations,
  that makes 36 primary specs: 32 ROBUST, 3 FRAGILE, 1 NOT SUPPORTED.
- **Abundance (11):** the agent found no change ("Nothing passes FDR < 0.05"). Tested type by
  type as "proportion differs", 4 are FRAGILE and 7 NOT SUPPORTED. After BH over the 11 types
  within each run *(ours)*, no type passes in more than 5 of 21 runs, and 17 would be needed to
  contradict the agent. This is consistent with the null, not evidence for it. On the agent's
  own labels, without reclustering, Recluster's test gives Megakaryocyte/platelet p = 0.0017,
  which passes BH over the 11 types (0.018 *(ours)*); the grid almost never reproduces that
  population.
- **Secondary (7):** 4 NOT SUPPORTED, 1 FRAGILE, 2 ROBUST.
- **The agent's summary disagrees with its own tables in four places.** For example, it names
  AMICA1 and NELL2 as down "in NK", but its NK table has no row for either.

Three cautions. First, the population results come from a re-run decided after the first
results had been seen (amendment 9). The first 11 population specs named no sample column, so
Scrublet scored all 24,562 cells as one batch and removed 4,826 of them at baseline instead of
349, including 51 of the 107 pDC cells. Re-run with the sample column and no other change to the
specs, pDC went from FRAGILE 16/21 to ROBUST 19/21; CD4 T memory stayed FRAGILE 14/21 and
Megakaryocyte/platelet NOT SUPPORTED 1/21. Both sets are reported below. Second, the min_genes
100 run is identical to the baseline in all 65 specs, so the baseline counts twice; without the
duplicate, two NOT SUPPORTED verdicts would be FRAGILE at 4/20 (see Limitations). Third, each
condition was captured on its own 10x lane, so a ROBUST stim-vs-ctrl result cannot separate
interferon from lane.

Columns marked *(ours)* are our own analysis of Recluster's per-run output (`provenance.json`),
made by `summarize.py` or read off the run records. They are not Recluster verdicts.

## Setup

### The run

- **Data.** Kang et al. 2018 (GEO GSE96583) as packaged by pertpy: 24,673 cells and 15,706
  genes from 8 donors, each with a control and an IFN-β-stimulated sample. The harness kept the
  raw counts, `label` (ctrl/stim) and `replicate` (donor), and removed the authors' cell types,
  clusters and embeddings.
- **Agent.** `claude-opus-5-5` through the Messages API, called by
  `recluster/bench/agent_audit.py`. No system prompt. One tool: a persistent Python 3.11 session
  in `run/workspace/` with the project's installed packages, `recluster` blocked from import and
  API keys removed from its environment. `max_tokens` 32,000, at most 60 API calls, 20 minutes per
  code cell. The prompt is [run/prompt.md](run/prompt.md) and the tool description is in
  [run/setup.json](run/setup.json). The prompt says the cells come from 8 donors and that
  `replicate` is the donor ID, because a user would describe their metadata. Neither file
  mentions pseudobulk, donors as the unit of analysis, doublets or time of day.
- **Size.** 16 API calls and 15 code cells, from 15:41:51 to 15:44:23 UTC on 12 November 2026
  (2 min 32 s, of which 61.8 s ran code). Tokens: 34 uncached input, 21,376 written to the cache,
  146,814 read from it, 9,523 output. The agent stopped on its own (`end_turn`)
  ([harness_log.txt](harness_log.txt), [run/run.json](run/run.json)).
- **No seed.** The API has no seed and this model does not accept `temperature`, so the run cannot
  be replayed exactly; a second run would give other clusters, labels and claims. The run was
  made once and did not fail, and it is kept unedited in [run/](run/). Its code can be re-executed.

### What the agent did

From [run/analysis.py](run/analysis.py) (all 15 cells, in order) and
[run/final_answer.md](run/final_answer.md):

- **QC.** Removed 111 cells with fewer than 200 genes (24,562 left) and genes detected in fewer
  than 10 cells (13,306 left). The matrix has no mitochondrial genes. No doublet detection.
- **Clustering.** Total-count normalisation, log1p, 2,000 highly variable genes, 30 principal
  components, Harmony on the condition (`label`), Leiden at resolution 0.8. Cells 005 to 008 patch
  a scipy 1.17 / scanpy incompatibility in the neighbour step.
- **Annotation.** From marker means per cluster: 12 clusters, two merged into "CD4 T naive", 11
  labels in `obs['cell_type']`. It wrote that two labels are less certain (CD4 T,
  Megakaryocyte/platelet).
- **Abundance.** Per-donor proportions, paired Wilcoxon signed-rank test per type, BH over the 11
  types.
- **Expression.** Pseudobulk per donor, condition and type (donors with at least 10 cells in both
  conditions, at least 3 donors), genes with at least 10 counts in `n_donors // 2 + 1` samples,
  PyDESeq2 with `~donor + condition`. pDC (1 usable donor) and Megakaryocyte/platelet (0) were
  skipped. Rows without an adjusted p were dropped. Table:
  `run/workspace/pseudobulk_DE_stim_vs_ctrl.csv` (9 types, 24,446 rows).
- It saved `annotated.h5ad` with the raw counts in `layers['counts']`, so the audits ran on that
  file directly (the fallback in protocol rule 5 was not needed).

### From the answer to claim specs

The protocol, [PROTOCOL.md](PROTOCOL.md), was written before the run. A statement about this
dataset becomes a spec when it maps onto a Recluster claim type: a gene up or down in a type
(`differential_expression`), a type more or less abundant (`composition_shift`), a type that
exists (`population_exists`). Populations are the agent's own labels, with no markers added.
Comparisons use condition `label`, stim against ctrl. Every spec except the first population run
(`P-*`, see amendment 9) names sample `[replicate, label]`, and all use the default grid and
thresholds. Amendments, each with its reason and timing in PROTOCOL.md:

1. Genes one statement names for one type go into one spec with `require: all`, because one spec
   per gene would have meant about 190 specs and about 11 hours; the pseudobulk p-values are then
   BH-adjusted over the named genes, which is stricter than one gene per spec.
2. A statement about a group of types is tested in each member where the agent's own table has
   the gene at adjusted p < 0.05 in the stated direction.
3. The abundance null is tested as "proportion differs" (`different`) for each of the 11 types; it
   is contradicted only by a ROBUST shift.
4. Only the genes the agent names are tested (17 of the 70 it says rise everywhere).
5. `make_specs.py` checks that every quote it maps appears verbatim in the final answer;
   `run_audits.py` runs `recluster audit` on each spec.
6. Group members that the agent's table does not support are tested too, as *secondary* specs
   (literal readings), and the missed caveat sentence is added the same way. Primary specs are
   the agent's statements mapped under the other rules; secondary specs are reported separately.
7. The abundance null is also judged with the agent's multiplicity: BH over the 11 types within
   each grid run; it is contradicted only if a type passes in at least 80 % of runs.
8. Each gene's own survival is reported next to the statement's verdict; the NK statement is split
   into one spec per gene; hedged labels are marked.
9. The 11 population specs (`P-*`) named only the label column, unlike the comparison specs, so
   doublet detection scored all cells as one batch. They are re-run as `P9-*` with sample
   `[replicate, label]` and no other change to the specs (the code and the worker count also
   changed; see below). `P9-*` replace `P-*` as the primary population results; both are
   reported. The amendment also corrects the wording of amendment 6: the
   agent's table has no row for those pairs, which does not show that the agent never tested
   them.

Amendments 1 to 5 were written after reading the answer and before the audits; the only
Recluster result seen by then was a timing check (ISG15 up in CD14 Mono, ROBUST 21/21, re-run in
the full set). Amendments 6 to 8 followed three extractions of the claims, made without seeing the
specs, and an adversarial review of `make_specs.py`, also before any audit. Amendment 9 was
written after the 54 audits had run and their results had been read, and before the re-run it
describes.

Every spec ran on the same 21-run grid: the baseline (min_genes 200, doublet removal on,
log-total normalisation, 2,000 HVGs, 30 PCs, resolution 1.0, seed 0), 14 one-at-a-time changes and
6 random combinations. Run numbers below are the `index` in `provenance.json`: 0 is the baseline,
1 to 14 the one-at-a-time changes, 15 to 20 the random combinations. In the five runs set to
min_genes 500, Recluster capped the threshold at 442 genes, which keeps QC cell loss under its
25 % limit; the run labels record this, and "min_genes 500" below means that capped setting. The
mitochondrial axis is absent (no mitochondrial genes) and so is the batch-correction axis (the
specs name no batch column): no grid run uses Harmony, which the agent did.

In each run a cluster stands for the agent's population when at least 50 % of its cells carry the
label. A population claim holds with at least 20 matched cells and Jaccard at least 0.50 to the
agent's cells; the Jaccard is computed on the agent's cells that survive the run's QC and doublet
removal. A composition claim holds when the paired t-test on per-donor logit proportions has
p < 0.05. An expression claim holds when every named gene has pseudobulk DESeq2 adjusted
p < 0.05 in the stated direction. The design is `~pair + condition` when every donor has a
pseudobulk in both conditions, and the unpaired `~condition` otherwise (below, under the cDC
claim, this matters). Recluster also tests each comparison claim once on the agent's labels as
given, with no reclustering (`as_labelled` in `provenance.json`); that result is recorded but is
not part of the verdict.

Every `provenance.json` records Recluster 0.6.0; it records the package version, not a commit.
The first 54 specs ran from 15:57 to 23:35 UTC on 12 November 2026 with 14 workers and took
27,487 s (7 h 38 min). Throughout, the checkout they imported was at 41d27d9, a local commit
with the same tree (the same files) as 9092e6f: the code of Recluster 0.6.0 plus
`recluster demo` and `--jobs`. The commits made between 22:17 and 23:20 UTC that day were made
in a separate git worktree; the checkout moved off 41d27d9 at 23:36 UTC, after the last spec had
finished, so all 54 ran on unchanged code. The 11 `P9-*` specs ran from 00:04 to 02:45 UTC on
13 November with 8 workers and took 9,684 s (2 h 41 min), on the package code of commit 1ffd61f.
In `recluster/` that code differs from 9092e6f by removed dead code and comments (CHANGELOG.md
records that the demo's key numbers hash identically before and after), the wording of the
collection-time check in reports, the demo's exit status, and running the sequential grid path
on one BLAS thread; the audits use worker processes (`--jobs`), which already ran on one BLAS
thread. Neither the code nor the worker count changed the grid runs here *(ours)*. In all 231
`P9-*` runs, the cells kept, doublets removed, number of clusters, matched clusters and matched
cell count equal those of the abundance spec for the same label, which names the same sample
column and ran on 41d27d9 with 14 workers. In the three runs without doublet removal, which the
sample column does not affect, every `P9-*` spec matches its `P-*` run exactly (same matched
cells and Jaccard).

## Results

### Populations (11 primary specs, re-run under amendment 9)

These are the `P9-*` specs, with sample `[replicate, label]`. *Median Jaccard (ours)* is over all
21 runs; a run in which no cluster matched counts as 0. *Removed as doublets at baseline (ours)*
is the number of the agent's cells for that label that the baseline run removed before matching
(QC removed none at baseline).

| Agent's label | Cells | Agent's hedge | Verdict | Runs held | Median Jaccard to the agent's cells *(ours)* | Removed as doublets at baseline *(ours)* |
|---|---|---|---|---|---|---|
| CD4 T naive | 6,650 | "the CD4 gene itself is barely detected" | ROBUST | 21/21 | 0.73 | 15 |
| CD14 Mono | 5,692 | | ROBUST | 21/21 | 0.96 | 143 |
| CD4 T memory | 3,429 | "the CD4 gene itself is barely detected" | **FRAGILE** | 14/21 | 0.57 | 18 |
| B | 2,696 | | ROBUST | 21/21 | 0.97 | 50 |
| NK | 1,756 | | ROBUST | 21/21 | 0.90 | 48 |
| CD8 T | 1,742 | | ROBUST | 20/21 | 0.76 | 32 |
| CD16 Mono | 1,093 | | ROBUST | 19/21 | 0.88 | 22 |
| T (stress/HSP-high) | 866 | | ROBUST | 21/21 | 0.73 | 4 |
| cDC | 436 | | ROBUST | 20/21 | 0.78 | 10 |
| pDC | 107 | | ROBUST | 19/21 | 0.96 | 7 |
| Megakaryocyte/platelet | 95 | "probably includes platelets stuck to T cells" | **NOT SUPPORTED** | 1/21 | 0.00 | 0 |

9 ROBUST, 1 FRAGILE, 1 NOT SUPPORTED. At baseline doublet detection removed 349 cells in all,
the same number as in the comparison specs, and at most 6.5 % of any label (7 of the 107 pDC
cells).

**The first run (`P-*`, superseded).** The first 11 population specs named only the label
column. Without a sample column, Scrublet scored all 24,562 cells as one capture and its
fallback removed 4,826 at baseline. That took between 7.4 % (7 of the 95 Megakaryocyte/platelet
cells) and 47.7 % (51 of the 107 pDC cells) of each label before matching: 46.3 % of cDC, 29.9 %
of CD4 T memory and 22.7 % of CD8 T, for example. Their verdicts, as they ran:

| Agent's label | `P-*` verdict | Runs held | `P9-*` verdict | Runs held |
|---|---|---|---|---|
| CD4 T naive | ROBUST | 21/21 | ROBUST | 21/21 |
| CD14 Mono | ROBUST | 21/21 | ROBUST | 21/21 |
| CD4 T memory | FRAGILE | 14/21 | FRAGILE | 14/21 |
| B | ROBUST | 21/21 | ROBUST | 21/21 |
| NK | ROBUST | 21/21 | ROBUST | 21/21 |
| CD8 T | ROBUST | 19/21 | ROBUST | 20/21 |
| CD16 Mono | ROBUST | 21/21 | ROBUST | 19/21 |
| T (stress/HSP-high) | ROBUST | 21/21 | ROBUST | 21/21 |
| cDC | ROBUST | 19/21 | ROBUST | 20/21 |
| pDC | **FRAGILE** | 16/21 | ROBUST | 19/21 |
| Megakaryocyte/platelet | NOT SUPPORTED | 1/21 | NOT SUPPORTED | 1/21 |

One verdict changed: pDC, from FRAGILE to ROBUST. In the first run no pDC cluster matched at
resolution 0.3 or 0.5, or in three of the five runs with min_genes 500 (one one-at-a-time, two
random). With the sample column it matched in all one-at-a-time runs, at Jaccard 0.79 to 0.99,
and failed only in random combinations 15 and 18 (both min_genes 500), where no cluster matched
*(ours)*. CD8 T, CD16 Mono and cDC moved by one or two runs and stayed ROBUST: CD16 Mono now
fails at resolution 0.3 (Jaccard 0.499) and in random combination 15 (0.42). Amendment 9 was
written after the `P-*` results had been read, so the pDC change should be weighed with that in
mind.

### Abundance (11 primary specs: the agent's null)

The agent: "Nothing passes FDR < 0.05." Each spec tests the opposite, that the type's proportion
differs between stim and ctrl, so a run that holds is a run in which that type, taken alone,
shifts at p < 0.05. The agent used a paired Wilcoxon signed-rank test on proportions; Recluster
uses a paired t-test on logit proportions. The agent's figures are from its
`abundance_stim_vs_ctrl.csv`.

| Type | Agent: ctrl → stim (mean %) | Agent's adjusted p | Recluster verdict on "differs" | Runs held | Not testable |
|---|---|---|---|---|---|
| CD14 Mono | 25.81 → 23.96 | 0.086 | FRAGILE | 15/21 | 0 |
| CD4 T naive | 25.37 → 24.67 | 0.634 | FRAGILE | 12/21 | 0 |
| CD4 T memory | 14.91 → 15.19 | 0.928 | FRAGILE | 5/21 | 4 |
| B | 10.09 → 10.51 | 0.393 | FRAGILE | 5/21 | 0 |
| NK | 6.94 → 7.03 | 0.945 | NOT SUPPORTED | 1/21 | 0 |
| CD8 T | 6.23 → 6.80 | 0.241 | NOT SUPPORTED | 4/21 | 1 |
| CD16 Mono | 4.38 → 4.80 | 0.393 | NOT SUPPORTED | 2/21 | 2 |
| T (stress/HSP-high) | 3.54 → 4.06 | 0.907 | NOT SUPPORTED | 1/21 | 0 |
| cDC | 2.14 → 1.76 | 0.241 | NOT SUPPORTED | 4/21 | 1 |
| pDC | 0.30 → 0.66 | 0.241 | NOT SUPPORTED | 0/21 | 2 |
| Megakaryocyte/platelet | 0.29 → 0.56 | 0.086 | NOT SUPPORTED | 0/21 | 20 |

Our per-run analysis of the same specs *(ours)*:

| Type | Held runs: lower / higher in stim | Runs passing BH over 11 types | Median per-run adjusted p |
|---|---|---|---|
| CD14 Mono | 15 / 0 | 5/21 | 0.23 |
| CD4 T naive | 9 / 3 | 5/21 | 0.16 |
| CD4 T memory | 2 / 3 | 5/21 | 0.32 |
| B | 0 / 5 | 1/21 | 0.50 |
| NK | 0 / 1 | 0/21 | 0.84 |
| CD8 T | 1 / 3 | 2/21 | 0.56 |
| CD16 Mono | 2 / 0 | 0/21 | 0.78 |
| T (stress/HSP-high) | 1 / 0 | 1/21 | 0.71 |
| cDC | 3 / 1 | 2/21 | 0.78 |
| pDC | 0 / 0 | 0/21 | 0.68 |
| Megakaryocyte/platelet | 0 / 0 | 0/21 | 1.00 |

Recluster: 4 FRAGILE, 7 NOT SUPPORTED, no ROBUST, so under amendment 3 the agent's null is not
contradicted for any type. Our per-run BH (amendment 7; a run that cannot test a type counts as
p = 1): no type passes in more than 5 of 21 runs, against the 17 needed, so the null is not
contradicted under the agent's own multiplicity either. NOT SUPPORTED here is consistent with the
null, not evidence for it: 8 donor pairs can only detect large shifts, which the agent also says.
With 8 pairs, the smallest two-sided p its exact Wilcoxon test can give is 0.0078; its
Megakaryocyte/platelet p is that value (0.0078125 in its CSV), and BH over 11 types makes it 0.086.
For pDC, the pooled-cell chi-square would call a shift in 19 of 21 runs and the donor-level test in
none (Recluster raises its pseudoreplication flag); the agent's donor-level test did not call it
either. On the agent's own labels the result differs for Megakaryocyte/platelet; see "Recluster on
the agent's own labels" below.

### Expression (25 primary specs)

*Per-gene runs held (ours)* counts, for each named gene, the runs in which that gene held at the
replicate level inside its spec (with the BH adjustment over the spec's genes), from the same
runs. A statement holds only when all its genes do, so each gene holds in at least as many runs
as its statement. "Each gene N" means no gene did better than the statement; in those rows every
failed run was one that could not test the claim.

> 70 genes go up significantly in all 9 tested types. Examples are ISG15, IFI6, IFIT1/2/3, MX1/MX2, IRF7, OAS1, LY6E, IFITM1/2/3, RSAD2, CMPK2, HERC5 and BST2.

17 genes up, in each of the 9 types the agent tested.

| Type | Verdict | Runs held | Not testable | Per-gene runs held *(ours)* |
|---|---|---|---|---|
| B | ROBUST | 21/21 | 0 | each gene 21 |
| CD14 Mono | ROBUST | 21/21 | 0 | each gene 21 |
| CD16 Mono | ROBUST | 19/21 | 2 | each gene 19 |
| CD4 T naive | ROBUST | 21/21 | 0 | each gene 21 |
| CD4 T memory | ROBUST | 17/21 | 4 | each gene 17 |
| CD8 T | ROBUST | 20/21 | 1 | each gene 20 |
| NK | ROBUST | 21/21 | 0 | each gene 21 |
| T (stress/HSP-high) | ROBUST | 21/21 | 0 | each gene 21 |
| cDC | ROBUST | 20/21 | 1 | each gene 20 |

> Examples are CCL8 (log2 fold change 9.8 in CD14 monocytes), CXCL10 (9.7 in CD14 monocytes, 12.0 in cDC), CXCL11, TNFSF10 and IFI27 (7.5), plus IL27, IDO1 and CH25H in CD16 monocytes.

The sentence follows "Monocytes and DCs respond most". The agent's table supports all five
chemokine and interferon genes in all three myeloid types, so there is no secondary spec.

| Genes | Type | Verdict | Runs held | Not testable | Per-gene runs held *(ours)* |
|---|---|---|---|---|---|
| CCL8, CXCL10, CXCL11, TNFSF10, IFI27 up | CD14 Mono | ROBUST | 21/21 | 0 | each gene 21 |
| CCL8, CXCL10, CXCL11, TNFSF10, IFI27 up | CD16 Mono | ROBUST | 19/21 | 2 | each gene 19 |
| CCL8, CXCL10, CXCL11, TNFSF10, IFI27 up | cDC | ROBUST | 20/21 | 1 | each gene 20 |
| IL27, IDO1, CH25H up | CD16 Mono | ROBUST | 19/21 | 2 | each gene 19 |

> Monocytes and cDC lose classical-monocyte and inflammatory genes: CD14 (−3.1 in CD16 monocytes), VCAN (−4.4), CLEC5A (−4.8), MMP9 (−4.3), IL8 (−4.1), CD93 and S100A8.

Primary: the genes the agent's table has down in each type. The rest are in the secondary table.

| Genes | Type | Verdict | Runs held | Not testable | Per-gene runs held *(ours)* |
|---|---|---|---|---|---|
| CD14, VCAN, CLEC5A, MMP9, IL8, CD93, S100A8 down | CD14 Mono | ROBUST | 21/21 | 0 | each gene 21 |
| CD14, IL8, S100A8 down | CD16 Mono | ROBUST | 19/21 | 2 | each gene 19 |
| CD14, VCAN, CLEC5A, IL8, S100A8 down | cDC | **FRAGILE** | 5/21 | 1 | S100A8 5, CD14 6, VCAN 17, CLEC5A 18, IL8 20 |

> Lymphocytes lose only a few genes, and modestly, e.g. AMICA1, NELL2, ITM2C and S100A4 in NK.

Tested in NK, one spec per gene (amendment 8), whatever the agent's table says.

| Gene | Type | Verdict | Runs held | Not testable | In the agent's NK table |
|---|---|---|---|---|---|
| AMICA1 down | NK | ROBUST | 19/21 | 0 | no row |
| NELL2 down | NK | **FRAGILE** | 6/21 | 0 | no row |
| ITM2C down | NK | ROBUST | 21/21 | 0 | log2FC −1.15, adjusted p 0.0022 |
| S100A4 down | NK | ROBUST | 21/21 | 0 | log2FC −1.85, adjusted p 5.8e-08 |

> CXCL10, CCL8 and CXCL11 also appear strongly "up" in T, B and NK cells.

The agent offers this as a caveat (it suspects RNA from monocytes), but it is still a statement
that the genes go up in these types. Primary: the gene-type pairs its table supports.

| Genes | Type | Verdict | Runs held | Not testable | Per-gene runs held *(ours)* |
|---|---|---|---|---|---|
| CXCL10, CCL8 up | CD4 T naive | ROBUST | 21/21 | 0 | each gene 21 |
| CXCL10, CCL8 up | CD4 T memory | ROBUST | 17/21 | 4 | each gene 17 |
| CXCL10 up | T (stress/HSP-high) | ROBUST | 21/21 | 0 | 21 |
| CXCL10, CCL8, CXCL11 up | B | ROBUST | 21/21 | 0 | each gene 21 |
| CXCL10, CCL8 up | NK | ROBUST | 21/21 | 0 | each gene 21 |

All 25 expression specs: 23 ROBUST, 2 FRAGILE.

### Secondary specs: literal readings (7 specs, amendment 6)

These test the members of a group statement that the agent's own table does not support (the
gene has no row for that type). They are not counted in the primary results.

| Statement | Genes | Type | Verdict | Runs held | Not testable | Per-gene runs held *(ours)* |
|---|---|---|---|---|---|---|
| Monocytes and cDC lose ... | VCAN, CLEC5A, MMP9, CD93 down | CD16 Mono | NOT SUPPORTED | 0/21 | 2 | VCAN 19, CLEC5A 19, MMP9 0, CD93 0 |
| Monocytes and cDC lose ... | MMP9, CD93 down | cDC | NOT SUPPORTED | 0/21 | 1 | MMP9 3, CD93 0 |
| CXCL10, CCL8 and CXCL11 ... | CXCL11 up | CD4 T naive | ROBUST | 21/21 | 0 | 21 |
| CXCL10, CCL8 and CXCL11 ... | CXCL11 up | CD4 T memory | FRAGILE | 15/21 | 4 | 15 |
| CXCL10, CCL8 and CXCL11 ... | CXCL10, CCL8, CXCL11 up | CD8 T | NOT SUPPORTED | 4/21 | 1 | CXCL10 20, CCL8 20, CXCL11 4 |
| CXCL10, CCL8 and CXCL11 ... | CCL8, CXCL11 up | T (stress/HSP-high) | NOT SUPPORTED | 1/21 | 0 | CCL8 21, CXCL11 1 |
| CXCL10, CCL8 and CXCL11 ... | CXCL11 up | NK | ROBUST | 20/21 | 0 | 20 |

4 NOT SUPPORTED, 1 FRAGILE, 2 ROBUST. The four failures come from three genes: MMP9 and CD93 in
CD16 Mono and cDC, and CXCL11 in CD8 T and stress-high T. The agent's table has no row for any of
these pairs. CXCL11 in CD4 T memory is FRAGILE, 15/21, with 4 runs that could not test it. Within
their specs, VCAN and CLEC5A in CD16 Mono, CXCL10 and CCL8 in CD8 T, and CCL8 in stress-high T hold
in 19 to 21 runs *(ours)*, although the agent's table has no row for them either. Recluster raises
the pseudoreplication flag on the CD8 T and stress-high T specs: the cell-level test supports them
in 20 and 21 of 21 runs.

### Recluster on the agent's own labels

For each comparison claim Recluster also runs the test once on the agent's labels as given, with
no reclustering and no QC or doublet removal beyond the agent's own (`as_labelled` in
`provenance.json`). It is not part of the verdict. The population specs have no such test.

- **Expression.** All 25 primary specs hold on the agent's labels. The cDC down-genes hold with
  the paired design, 8 donors against 8; the largest adjusted p is 2.2e-4 (VCAN). In the agent's cDC
  cells CD14 is detected in 22.6 % of control and 6.4 % of stimulated cells. NELL2 in NK holds
  (log2FC −2.20, adjusted p 0.029), and so does AMICA1 (−2.26, 0.0041).
- **Secondary.** The four NOT SUPPORTED specs also fail on the agent's labels. CXCL11 has adjusted
  p 0.61 in CD8 T and 0.23 in stress-high T; MMP9 and CD93 have 0.91 in CD16 Mono, and 0.14 and
  0.88 in cDC. The other three secondary specs hold.
- **Abundance.** Three types pass p < 0.05 on their own: Megakaryocyte/platelet (p = 0.0017,
  higher in stim), CD14 Mono (0.030, lower) and CD8 T (0.049, higher). BH over the 11 types
  *(ours)* gives 0.018, 0.16 and 0.18, so Megakaryocyte/platelet passes. On the agent's own
  labels, then, Recluster's test with the agent's multiplicity would reject the null for that
  type. The test rests on 2 to 20 Megakaryocyte/platelet cells per sample. The grid does not
  reject the null because the population is almost never reproduced: 20 of 21 runs cannot test
  it and count as p = 1.

## What broke the claims that were not ROBUST

Read from each claim's `breakers` in `out/results.json` and the run records in its
`provenance.json`. The per-claim output folders (`provenance.json`, `runs.csv`, `report.html`)
are not in the repository; `out/results.json`, `out/results_tables.md` and `out/summary.json`
are. Four of the reports are published and linked below. The population claims are the `P9-*`
specs.

**Megakaryocyte/platelet exists: NOT SUPPORTED, 1/21**
([report](https://pranakkin.github.io/Recluster/reports/agent/agent-megakaryocyte-platelet-is-a-distinct-population-per-sa.html)).
In 20 of 21 runs no cluster had at least half of its cells labelled Megakaryocyte/platelet, so
nothing matched; that includes the three runs with doublet removal off. It held only at
resolution 0.3 (157 cells, Jaccard 0.57). It also fails at baseline, where doublet detection
removed none of its 95 cells, and under 13 of the 14 one-at-a-time changes. The first run (`P-*`)
gave the same result: 1/21, held only at resolution 0.3 (143 cells, Jaccard 0.59).
The agent hedged this label: "this cluster also shows CD3D, so it probably includes platelets
stuck to T cells." Its own marker means for the cluster (cell 009) show PPBP 4.47 and PF4 3.07
next to CCR7 1.81, SELL 1.35, IL7R 0.85 and CD3E 0.59. It also wrote "This rests on only 31 vs 64
cells" about the abundance trend. The hedge is not carried into the main findings, which list
"more platelets (0.29% → 0.56%)" as one of the two strongest trends.

**CD4 T memory exists: FRAGILE, 14/21.** It holds at baseline (Jaccard 0.58; doublet detection
removed 18 of its 3,429 cells) and fails in 7 runs. Among the one-at-a-time changes, no cluster
matched at resolution 0.3 or 0.5; with 50 PCs and with seed 1 the matched clusters held only 920
and 929 cells (Jaccard 0.19 in both); with doublet removal off the matched set reached Jaccard
0.46, under 0.50. Two random combinations with min_genes 500 failed: run 15 (which includes
resolution 0.3) with no cluster matched, and run 18 (with 50 PCs) at Jaccard 0.28. With 1,000
HVGs it held at Jaccard 0.501, just over the line. The median Jaccard over all runs is 0.57
*(ours)*: the agent's CD4 T memory cells overlap only partly with the clusters matched to them
on the grid. The first run (`P-*`) also held in 14 of 21, but not in the same runs: there the
1,000-HVG run failed (Jaccard 0.46) and the seed 1 run held. The agent hedged both CD4 T labels
because "the CD4 gene itself is barely detected"; CD4 T naive is ROBUST.

**CD14, VCAN, CLEC5A, IL8, S100A8 down in cDC: FRAGILE, 5/21**
([report](https://pranakkin.github.io/Recluster/reports/agent/agent-myeloid-down-in-cdc.html)).
It fails at baseline: CD14 has pseudobulk adjusted p 0.052 and S100A8 0.43. It holds with
Pearson-residual normalisation, with 50 PCs and in three random combinations, and fails under
every other one-at-a-time change; at resolution 0.3 no cDC cluster matched. VCAN, CLEC5A and IL8
hold in 17, 18 and 20 runs within the spec; CD14 in 6 and S100A8 in 5 *(ours)*. In the five runs where the statement held, the matched control cDC set
had 218 to 317 cells and CD14 was detected in 20.4 % to 38.5 % of them; in the 15 other testable
runs it had 83 to 201 cells and 1.6 % to 9.1 % *(ours)*. On the agent's own cDC cells the figure
is 22.6 %, and the statement holds there. This suggests the decrease depends on where the cDC
boundary falls, and in particular on CD14-expressing cells being in the control cDC set. We did
not test this further. The design differs too: 15 of the 21 cDC runs, the baseline among them,
used the unpaired `~condition`, because not every donor had a cDC pseudobulk in both conditions
(at baseline patient_1039 had none in stim). The agent used `~donor + condition` over the 7
donors with at least 10 cells in both. The design does not explain the failures on its own: 4 of
the 5 runs that held were unpaired (runs 4, 8, 16 and 19), and 4 runs that failed were paired
(runs 2, 3, 15 and 18). The agent's own cDC table has all five genes down: CD14 log2FC −3.21
(adjusted p 6.6e-08), VCAN −4.44 (1.0e-07), CLEC5A −4.59 (5.8e-08), IL8 −2.67 (3.6e-08), S100A8
−3.55 (1.6e-08), from 7 donors.

**NELL2 down in NK: FRAGILE, 6/21**
([report](https://pranakkin.github.io/Recluster/reports/agent/agent-nk-down-in-nk-nell2.html)).
At baseline the pseudobulk log2FC is −2.22 with adjusted p 0.062, just above 0.05; NELL2 is detected
in 0.11 % of stimulated and 1.33 % of control NK cells, so it is close to absent in NK and the
verdict says little either way. It holds with min_genes 500, with doublet removal off, with 1,000
HVGs and in three random combinations (runs 15, 16 and 20), and fails in the other 15 runs, the
baseline among them. On the agent's own labels it holds (adjusted p 0.029). The cell-level test
supports it in 21 of 21 runs, so Recluster raises the pseudoreplication flag. The agent's NK table
has no NELL2 row; NELL2 is in its CD4 T naive table (log2FC −1.69, adjusted p 2.1e-22) and CD4 T
memory table (−1.39, 5.2e-04) only.

**Abundance, the four FRAGILE types.** Here a run that holds is a nominal shift against the
agent's null, so "broke" means the shift disappeared.

- CD14 Mono, 15/21
  ([report](https://pranakkin.github.io/Recluster/reports/agent/agent-null-cd14-mono-proportion-differs-stim-vs-ctrl.html)):
  lower in stim at baseline (paired t-test p = 0.029) and in all 15 runs that held *(ours)*. It
  fails in 6 runs: 4 of the 5 with Pearson residuals (the one-at-a-time run and three random
  combinations), resolution 1.5 (p = 0.051), and one random combination with resolution 1.5 and
  doublet removal off. Only 5 runs pass BH over the 11 types *(ours)*. The agent found the same
  trend (lower in 7 of 8 donors, p = 0.016, adjusted p = 0.086). Recluster's benchmark gives the
  same claim, on the authors' labels, FRAGILE 11/21 ([bench/RESULTS.md](../RESULTS.md), case 6).
  It scores that result as a false alarm, because the case was registered as having no known
  shift.
- CD4 T naive, 12/21: held at baseline (p = 0.044), but 3 of the 12 runs that held go the other
  way *(ours)*. Lost at resolution 0.3, 0.5 and 1.5, with 3,000 HVGs, with 50 PCs, with seed 1 or
  2, and in two random combinations (runs 15 and 17).
- CD4 T memory, 5/21: not at baseline (p = 0.077); held with doublet removal off, Pearson
  residuals, 1,000 HVGs and in two random combinations, in both directions *(ours)*. 4 runs could
  not test it: in three no cluster matched, in one the matched set was in only 1 control sample.
- B, 5/21: held only in the five runs with min_genes 500, all higher in stim (log-odds 0.16 to
  0.17) *(ours)*.

The seven NOT SUPPORTED types are consistent with the agent's null. CD8 T holds at baseline
(p = 0.0038, higher in stim), in the min_genes 100 run that duplicates it, and with 3,000 HVGs
and with 20 PCs. The Megakaryocyte/platelet composition spec could not test 20 of its 21 runs.

## Where the agent's summary disagrees with its own tables

Each checked against `run/workspace/*.csv`, the printed cell outputs in `run/cells/` and
`run/final_answer.md`.

1. **AMICA1 and NELL2 "in NK".** The NK table (1,948 genes) has no row for either. AMICA1 is down
   in CD4 T naive (−1.74), CD4 T memory (−1.76), CD8 T (−1.78), CD14 Mono (−1.50) and cDC (−4.36);
   NELL2 in CD4 T naive (−1.69) and CD4 T memory (−1.39). The NK down-list the agent printed (cell
   013) is RAMP1, S100A4, TECR, ITM2A, SRM, IMPDH2, FGFBP2, CXXC5: of the four genes it named,
   only S100A4 is there. ITM2C is in its NK table (−1.15, adjusted p 0.0022) but appears in the
   printed lists of CD4 T naive and CD8 T. Under Recluster, AMICA1 in NK is ROBUST (19/21) and
   NELL2 FRAGILE (6/21).
2. **"Everything else, including B, NK, T and DC proportions, changes by less than about 0.5
   percentage points."** Its abundance table has CD4 T naive 25.37 % → 24.67 %, a change of −0.70
   points, which is clearly over. CD8 T 6.23 % → 6.80 % (+0.57) and stress-high T 3.54 % → 4.06 %
   (+0.52) are close to "about 0.5". The other six are within 0.42 points.
3. **"CXCL10, CCL8 and CXCL11 also appear strongly "up" in T, B and NK cells."** Of the 18
   gene-type pairs (3 genes, 6 lymphocyte types), 8 have no row in its table: CXCL10 in CD8 T;
   CCL8 in CD8 T and stress-high T; CXCL11 in CD4 T naive, CD4 T memory, CD8 T, stress-high T and
   NK. CXCL11 has a row in B only. The table the agent printed in cell 013 shows NaN for CXCL10 in
   CD8 T and for CCL8 in CD8 T and stress-high T. The 10 pairs with rows are all up at adjusted
   p < 0.05.
4. **"1,364 changed genes in CD14 monocytes versus 141–324 in lymphocytes."** Its own count table
   gives stress-high T 118 up and 12 down, 130 in all, below 141.

## Statements that were not given specs

Quotes are exact. No Recluster claim type covers these, unless a reason is given:

- The two abundance trends, "lower in 7 of 8 donors" (CD14 Mono) and "higher in 8 of 8 donors"
  (Megakaryocyte/platelet). The agent presented them as not significant, so they were not given
  directional specs; the "differs" specs cover both types. All 15 CD14 Mono runs that held were
  lower in stim *(ours)*, and Megakaryocyte/platelet held in none.
- QC counts (111 cells removed, 24,562 left), the absence of mitochondrial genes, and that no
  doublet detection was run: statements about its own procedure.
- The choice of Harmony on the condition, "Without that correction, interferon dominates the data
  and the same cell type ends up in separate stim and ctrl clusters", "Every cluster contains both
  conditions in similar numbers" and "Each cell type has a matching stim and ctrl cluster after
  correction": a method choice and properties of its clustering, not of the data.
- The marker table used for annotation: the population specs use the agent's labels only.
- "pDC and Megakaryocyte/platelet were too small to test": a statement about its own filters.
- DE gene counts per type (for example 812 up and 552 down in CD14 Mono), "70 genes" up in all 9
  types, "581 up-regulated genes" found only in monocytes or cDC, and 1,364 changed genes in CD14
  Mono: counts from its own thresholds.
- "a very large gene-expression response in every cell type I tested": the size is not testable;
  the direction is covered by the shared interferon-gene specs.
- Fold-change magnitudes (ISG15 about 25× to 200×, CCL8 log2FC 9.8, CXCL10 12.0 in cDC, and the
  others): Recluster tests direction, not size.
- Detection percentages (ISG15 "detected in 97–100% of stimulated cells", CCL8 "in only about 11%
  of stimulated T cells versus 90% of CD14 monocytes").
- Comparisons between types ("Monocytes and DCs respond most", "The lymphocyte response is smaller
  and mostly increases").
- The explanation of the lymphocyte chemokine signal as "RNA leaking from monocytes or
  monocyte–lymphocyte doublets": a mechanism.
- The power statement ("With only 8 donor pairs, the test can only detect fairly large shifts").
- The scipy patch and its claim that it "doesn't change the results".

## Limitations

- **One run, one dataset, one prompt.** A second run would cluster, label and claim differently;
  the counts describe this run only.
- **Amendments after reading the answer.** All nine amendments were written after the agent's
  answer was read; PROTOCOL.md says for each when it was written relative to the run and the
  audits and why. Amendments 2 and 6 decide which type-gene pairs are primary; amendment 7
  decides how the abundance null is judged. Each is a choice, and other choices would change the
  counts. Amendment 6 was added because amendment 2 favoured the agent.
- **A re-run after seeing results.** Amendment 9 was written after the 54 audits had been read.
  It fixes an error in `make_specs.py`: the population specs named no sample column, although
  the comparison specs did, so Scrublet scored all cells as one batch and removed 4,826 cells at
  baseline instead of 349. The fix changed one verdict, pDC, from FRAGILE 16/21 to ROBUST 19/21,
  in the agent's favour. Both runs are reported under Populations; the re-run counts as primary.
- **Grouped genes.** Statements are tested as groups with `require: all` to cut run time
  (amendment 1). One weak gene fails the whole statement (S100A8 and CD14 in cDC), and the BH
  adjustment over the named genes is stricter than single-gene specs.
- **One grid run duplicates the baseline.** The agent had already removed cells with fewer than
  200 genes, so min_genes 100 (run 1) removes nothing more and gives the same result as the
  baseline in all 65 specs; Recluster did not merge the two. Counting the baseline once *(ours,
  not a re-run)* would change two verdicts at the 20 % line: cDC abundance from NOT SUPPORTED to
  FRAGILE, 4/20 (still not a contradiction of the null under amendment 3), and the secondary
  CD8 T spec from NOT SUPPORTED to FRAGILE, 4/20. No other verdict in the 65 specs changes. The
  verdicts above are Recluster's, as run.
- **What the grid does not vary.** There is no batch column, so no grid run uses batch correction.
  The agent's Harmony on the condition is not among the choices tested; Harmony on the condition
  can remove part of the condition effect, and the audit says nothing about that choice. Recluster
  matches the agent's labels to its own clusters; a label that no cluster carries by majority
  fails, even if the cells exist.
- **ROBUST is not truth.** It means the claim survives the grid, tested across donors. In Kang 2018
  each condition was captured on its own 10x lane, so a lane artefact would also come out ROBUST.
- **Time of day.** The prompt never mentions collection time. Recluster's time check ran on the
  43 comparison specs (population specs have none) and returned NO TIME CONFOUND DETECTED for all
  of them: two clock genes (CRY1, DBP) differ between stim and ctrl, but not in the pattern a time
  shift produces.

## How to reproduce

1. **The agent run:** `uv run python -m recluster.bench.agent_audit` (needs `ANTHROPIC_API_KEY`;
   refuses to overwrite an existing `run/`). It will not reproduce this run: the API has no seed.
2. **Inputs for this run's audits.** The `.h5ad` files in `run/workspace/` are not committed
   (`.gitignore`). `pbmc_ifnb.h5ad` is rebuilt from `data/kang_2018.h5ad` by `prepare_input()` in
   `recluster/bench/agent_audit.py` (SHA-256 in `run/setup.json`); `annotated.h5ad` by re-executing
   `run/analysis.py` in `run/workspace/` (SHA-256 `fcb5f9b5...` in every `provenance.json`). We have
   not checked that re-execution reproduces that checksum.
3. `uv run python bench/agent_audit/make_specs.py`: writes `claims/*.yaml` and
   `claims/claims.json` (65 specs: the 11 `P9-*` population specs, 11 abundance and 25
   expression specs that are primary, 7 secondary, and the 11 superseded `P-*` specs; each
   `P9-*` entry names the spec it `supersedes`).
4. `uv run python bench/agent_audit/run_audits.py --jobs 14`: runs `recluster audit` on each spec
   into `out/` and writes `out/summary.json`; resumes where it stopped. Here the first 54 specs
   took 7 h 38 min with `--jobs 14` and the 11 `P9-*` specs 2 h 41 min with `--jobs 8` (macOS
   arm64).
5. `uv run python bench/agent_audit/summarize.py`: writes `out/results.json` and
   `out/results_tables.md`, with the per-gene counts and the per-run BH.

The audit logs are [audit_log.txt](audit_log.txt) (the first 54 specs) and
[audit_log_a9.txt](audit_log_a9.txt) (the 11 `P9-*` specs).
