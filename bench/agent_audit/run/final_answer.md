I annotated 11 cell types. Interferon-beta causes a very large gene-expression response in every cell type I tested, but **no cell type changes significantly in abundance**.

## What I did
- **Quality control:** I removed 111 cells with fewer than 200 genes, leaving 24,562 cells (8 donors, each with a control and a stimulated sample). The data contain no mitochondrial genes, so I couldn't filter on that, and I didn't run doublet detection.
- **Clustering:** I normalised and log-transformed the counts, used 2,000 variable genes and 30 principal components, and then corrected for condition with Harmony. Without that correction, interferon dominates the data and the same cell type ends up in separate stim and ctrl clusters. Every cluster contains both conditions in similar numbers.
- **Annotation:** I named the clusters from known marker genes.

| Cell type | Markers | Cells |
|---|---|---|
| CD4 T naive (2 clusters merged) | CD3D, CCR7, SELL, IL7R | 6,650 |
| CD14 monocytes | CD14, LYZ, S100A8 | 5,692 |
| CD4 T memory | CD3D, IL7R, CREM, some CCL5 | 3,429 |
| B | MS4A1, CD79A | 2,696 |
| NK | GNLY, NKG7, GZMB, no CD3 | 1,756 |
| CD8 T | CD8A/B, CCL5, GZMH | 1,742 |
| CD16 monocytes | FCGR3A, MS4A7 | 1,093 |
| T cells with high stress genes | CD3, IL7R, CCR7 plus heat-shock genes (DNAJB1, CACYBP, HSPA1A) | 866 |
| cDC | FCER1A, CLEC10A, HLA-DR high | 436 |
| pDC | LILRA4, IL3RA, GZMB | 107 |
| Megakaryocyte/platelet | PPBP, PF4 | 95 |

Two labels are less certain:
- **CD4 T cells:** the CD4 gene itself is barely detected, so I called these CD4 because they are CD3-positive and CD8-negative.
- **Megakaryocyte/platelet:** this cluster also shows CD3D, so it probably includes platelets stuck to T cells.

## Changes in abundance (stim vs ctrl)
I compared each donor's cell-type proportions between conditions (paired Wilcoxon test, 8 donors, with false-discovery-rate correction). **Nothing passes FDR < 0.05.** Two trends are worth noting:
- **CD14 monocytes:** 25.8% → 24.0%, lower in 7 of 8 donors (p = 0.016, adjusted p = 0.086).
- **Megakaryocyte/platelet:** 0.29% → 0.56%, higher in 8 of 8 donors (p = 0.008, adjusted p = 0.086). This rests on only 31 vs 64 cells.
- Everything else, including B, NK, T and DC proportions, changes by less than about 0.5 percentage points.

With only 8 donor pairs, the test can only detect fairly large shifts.

## Gene expression changes
I summed counts per donor, cell type and condition, then tested with DESeq2 while accounting for donor. Counts below are genes with adjusted p < 0.05 and at least a 2-fold change:

| Cell type | Up | Down |
|---|---|---|
| CD14 Mono | 812 | 552 |
| CD16 Mono | 369 | 135 |
| cDC | 332 | 141 |
| CD4 T naive | 247 | 77 |
| B | 244 | 74 |
| CD4 T memory | 211 | 21 |
| NK | 178 | 36 |
| CD8 T | 129 | 12 |
| Stress-high T | 118 | 12 |

pDC and Megakaryocyte/platelet were too small to test (fewer than 10 cells per sample in most donors).

- **A shared interferon response in every cell type:** 70 genes go up significantly in all 9 tested types. Examples are ISG15, IFI6, IFIT1/2/3, MX1/MX2, IRF7, OAS1, LY6E, IFITM1/2/3, RSAD2, CMPK2, HERC5 and BST2. ISG15 rises about 2⁴·⁵ to 2⁷·⁷ fold (roughly 25× to 200×), and is detected in 97–100% of stimulated cells.
- **Monocytes and DCs respond most:** they have the most changed genes and the largest fold changes. Examples are CCL8 (log2 fold change 9.8 in CD14 monocytes), CXCL10 (9.7 in CD14 monocytes, 12.0 in cDC), CXCL11, TNFSF10 and IFI27 (7.5), plus IL27, IDO1 and CH25H in CD16 monocytes. 581 up-regulated genes are found only in monocytes or cDCs.
- **Genes that go down are cell-type specific:**
  - Monocytes and cDC lose classical-monocyte and inflammatory genes: CD14 (−3.1 in CD16 monocytes), VCAN (−4.4), CLEC5A (−4.8), MMP9 (−4.3), IL8 (−4.1), CD93 and S100A8.
  - Lymphocytes lose only a few genes, and modestly, e.g. AMICA1, NELL2, ITM2C and S100A4 in NK.

**Caveat:** CXCL10, CCL8 and CXCL11 also appear strongly "up" in T, B and NK cells. CCL8, for example, is found in only about 11% of stimulated T cells versus 90% of CD14 monocytes. That suggests at least part of this signal is RNA leaking from monocytes or monocyte–lymphocyte doublets rather than the lymphocytes' own expression.

## Main findings
1. **Interferon-beta changes gene expression, not cell-type proportions.** Each cell type has a matching stim and ctrl cluster after correction. No abundance change is significant; the strongest trends are fewer CD14 monocytes (25.8% → 24.0%) and more platelets (0.29% → 0.56%), both adjusted p = 0.086.
2. **A shared set of 70 interferon-response genes rises in all 9 tested cell types**, with ISG15 up about 25–200×.
3. **Monocytes and cDCs respond the most:** 1,364 changed genes in CD14 monocytes versus 141–324 in lymphocytes. Their up-regulated genes include chemokines (CCL8, CXCL10, CXCL11, log2 fold change about 10), and they lose CD14, VCAN, CLEC5A and IL8.
4. **The lymphocyte response is smaller and mostly increases**, e.g. CD8 T has 129 genes up and 12 down. The chemokine signals in lymphocytes are likely partly leaked monocyte RNA.

A small fix to the software environment (scipy 1.17 vs scanpy) was needed for the neighbour-graph step; it doesn't change the results.

Files are in the working directory:
- `annotated.h5ad`
- `pseudobulk_DE_stim_vs_ctrl.csv`
- `abundance_stim_vs_ctrl.csv`
- `umap.png`
