# Design measurements

Numbers cited in `bench/RESULTS.md` ("What changed between the runs") and `docs/HOW_IT_WORKS.md` (section 8), recomputed from `data/` by `scripts/design_measurements.py` with the package's own functions. Unrounded values and per-item details are in `design_measurements.json`. "Agrees" means the measured value rounds to the cited one at the cited precision (or lies inside a cited bound).

Recluster 0.5.0, run 2026-11-11 11:15, 1.7 min. The cluster-level doublet mixture analysis (RESULTS item 9) is not part of this script.

| Item | Quantity | Cited | Measured | Agrees | How computed |
|---|---|---|---|---|---|
| 1 (RESULTS 2) | Stim-monocyte cluster whole-dataset AUROC, baseline -> min_genes=500 uncapped | 0.83 -> 0.79 | 0.834 -> 0.788 (automatic Scrublet threshold only: 0.829 -> 0.788) | yes | Pipeline.run at min_genes 200 and 500 with GridConfig(qc_max_cell_loss=1.0); cluster with the most stim cells labelled CD14+ Monocytes; auroc_by_group of the CD14/LYZ/S100A8 marker_score over the run's cells |
| 1 (RESULTS 2) | Same cluster, stratified within-sample AUROC (the current statistic) | not cited | 0.938 -> 0.945 | – | within_sample_auroc on the same runs, sample = replicate x label |
| 1 (RESULTS 2) | Fraction of Kang cells min_genes=500 removes | 44 % | 44.1 % | yes | share of cells with fewer than 500 detected genes in the raw counts |
| 1 (RESULTS 2) | min_genes=500 after the default 25 % cap | 441 | 441 | yes | Pipeline.effective_qc(PipelineParams(min_genes=500)) with GridConfig() |
| 2 (RESULTS 2) | Rejection rate at true AUROC 0.85: minimum over samples vs stratified | 95 % vs 0 % | 95.0 % vs 0.0 % | yes | 16 samples x (12 cluster + 300 other cells), binormal scores with equal variance, 200 replicates (seed 0 at 0.85, seed 1 at 0.90); rejection = statistic < 0.80 |
| 2 (RESULTS 2) | Mean of the stratified statistic at true AUROC 0.85 | 0.851 | 0.851 | yes | within_sample_auroc, mean over the 200 replicates |
| 2 (RESULTS 2) | Rejection rate of the minimum at true AUROC 0.90, 200 replicates | not cited (200 replicates vary; see the next row) | 33.0 % (stratified: 0.0 %) | – | as above, minimum of auroc_two_sets over the 16 samples |
| 2 (RESULTS 2) | Rejection rate of the minimum at true AUROC 0.90, 2,000 replicates | about 32 % | 31.9 % | yes | as above, 2,000 replicates, seed 2 |
| 3 (RESULTS 4) | CD4 T matched in every Kang grid run, Jaccard range | every run, 0.90-0.97 | 21 of 21 runs, 0.90-0.97 | yes | bench/out/kang/isg15-up-in-cd4-t-cells-stim-vs-ctrl-kang-2018/provenance.json: runs[].match (candidate_rule=either; candidates qualified by both, labels) |
| 4 (RESULTS 6) | Testable cell types; SE shrink, unpaired SE / paired SE | 7 types; 2- to 9-fold | 7 types; 2.1- to 8.8-fold | yes | evaluate_composition per cell_type (all cells as the run), sample = [replicate, label], without and with pair_key; SE of the logit difference: Welch vs per-donor differences |
| 4 (RESULTS 6) | Unpaired Welch p over the testable types | all >= 0.67, none called | min 0.67, 0 called | yes | p_logit_ttest of the unpaired evaluate_composition |
| 4 (RESULTS 6) | Paired calls; CD14+ Monocytes p, share stim vs ctrl, donors lower in stim | CD14+ Monocytes only; p = 0.035, 23.8 % vs 25.4 %, 7 of 8 | CD14+ Monocytes; p = 0.035, 23.8 % vs 25.4 %, 7 of 8 | yes | paired evaluate_composition (paired t on per-donor logit differences) |
| 4 (RESULTS 6) | Monocyte DE: share of unpaired calls kept under ~pair + condition; sign flips | 95 %; 0 | 95.4 % (3286 of 3444); 0 flips | yes | pseudobulk_de on cells labelled CD14+ Monocytes, 10086 genes with >= 10 counts; call = DESeq2 padj < 0.05; kept = also called with the same sign |
| 5 (RESULTS 8) | Expected doublet rate vs demuxlet, per GEO lane (ambiguous excluded) | within 0.7 percentage points | ctrl 11.6 % vs 11.5 %; stim 11.5 % vs 11.5 %; largest gap 0.13 pp | yes | load_batch2 + Pipeline.qc_stage at baseline QC, Scrublet per lane; expected = 0.8 % per 1,000 scored cells of the lane; demuxlet rate over singlets + doublets |
| 5 (RESULTS 8) | Planted PBMC3k doublets: median Scrublet score of injected cells | 0.10 | 0.099 | yes | Pipeline.qc_stage at baseline QC on data/planted/pbmc3k_doublets.h5ad (one subset) |
| 5 (RESULTS 8) | Share of injected cells caught by an oracle cutoff at the true count | 28 % | 28.0 % (top 150 scores) | yes | injected share among the k highest scores, k = injected cells passing QC |
| 6 (RESULTS 12) | RPS26 donor-to-donor log2 CPM range, control cells (control monocytes) | about 2.6 (kang.py comment: 1.7 in monocytes) | 2.57 (1.74) | yes | log2(CPM + 1) of per-donor pseudobulks of the ctrl cells, max - min over 8 donors |
| 6 (RESULTS 12) | SQSTM1 donor-to-donor log2 CPM range, control monocytes (control cells) | about 0.4 (kang.py comment: 0.6 over all control cells) | 0.41 (0.57) | yes | as above |
| 6 (RESULTS 12) | Planted SQSTM1 at half-amplitude 0.4: unadjusted pseudobulk padj | 0.065 | 0.065 (log2FC 0.41; genome-wide padj 1.000); at 0.6: 0.0070 (log2FC 0.66) | yes | plant_time_signal(seed=0) with TIME_SIGNAL['SQSTM1'] half-amplitude replaced; pseudobulk_de on cells labelled CD14+ Monocytes, sample = replicate, B vs A (4 vs 4), ~condition; padj over the one claimed gene |
| 7 (HOW_IT_WORKS 8) | Clock-gene per-gene BH padj, DBP | 0.012 | 0.012 | yes | exchangeable_units + panel_check on all Kang cells, sample = [replicate, label], stim vs ctrl (8 vs 8 units) |
| 7 (HOW_IT_WORKS 8) | Clock-gene per-gene BH padj, CRY1 | 0.037 | 0.037 | yes | exchangeable_units + panel_check on all Kang cells, sample = [replicate, label], stim vs ctrl (8 vs 8 units) |
| 7 (HOW_IT_WORKS 8) | Clock-gene per-gene BH padj, PER1 | 0.051 | 0.051 | yes | exchangeable_units + panel_check on all Kang cells, sample = [replicate, label], stim vs ctrl (8 vs 8 units) |
| 7 (HOW_IT_WORKS 8) | Clock-state index permutation p | 0.38 | 0.381 (12870 exact permutations) | yes | exchangeable_units + panel_check on all Kang cells, sample = [replicate, label], stim vs ctrl (8 vs 8 units) |
