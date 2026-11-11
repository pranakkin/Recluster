| # | Case | Claim | Expected | Verdict | Held | Cell-level | Outcome |
|---|---|---|---|---|---|---|---|
| 1 | ISG15 up in CD14+ monocytes after IFN-beta | ISG15 up in CD14+ monocytes, stim vs ctrl (Kang 2018) | ROBUST | **FRAGILE** | 16/21 | 100% | MISSED |
| 2 | Four interferon-stimulated genes up in CD14+ monocytes | ISG15, IFI6, IFIT3 and MX1 all up in CD14+ monocytes, stim vs ctrl (Kang 2018) | ROBUST | **FRAGILE** | 16/21 | 100% | MISSED |
| 3 | ISG15 up in CD4 T cells after IFN-beta | ISG15 up in CD4 T cells, stim vs ctrl (Kang 2018) | ROBUST | **FRAGILE** | 15/21 | 81% | MISSED |
| 4 | CD14+ monocytes exist across donors | CD14+ monocytes exist (Kang 2018) | ROBUST | **FRAGILE** | 16/21 | 76% | MISSED |
| 5 | NK cells exist across donors | NK cells exist (Kang 2018) | ROBUST | **FRAGILE** | 15/21 | 71% | MISSED |
| 6 | CD14+ monocyte proportion differs stim vs ctrl | CD14+ monocytes differ in abundance, stim vs ctrl (Kang 2018) | NOT SUPPORTED (no known shift; reported as measured) | **FRAGILE** | 8/21 | 100% ⚑ | FALSE ALARM |
| 7 | False: ISG15 differs between fake donor groups (control cells only) | FALSE: ISG15 differs in CD14+ monocytes between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 1/21 | 81% ⚑ | correct |
| 8 | False: LYZ differs between fake donor groups (control cells only) | FALSE: LYZ differs in CD14+ monocytes between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 100% ⚑ | correct |
| 9 | False: CD14+ monocyte proportion differs between fake donor groups | FALSE: CD14+ monocyte proportion differs between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 19% | correct |
| 10 | False: ISG15 up in fake group A (condition shuffled at donor level) | FALSE: ISG15 up in CD14+ monocytes in fake group A vs B (Kang, donor-level shuffle) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 100% ⚑ | correct |
| 11 | Injected synthetic doublets forming a fake population | Planted doublets: a CD14+ LYZ+ MS4A1+ CD79A+ population exists (PBMC3k) | not ROBUST; holds only when doublet removal is off | **ROBUST** | 23/23 | 100% | MISSED |
| 12 | Over-clustering that creates a spurious population | Over-clustering: cluster_14 is a distinct population (PBMC3k) | not ROBUST | **NOT SUPPORTED** | 0/23 | 0% | caught |
| 13 | Control for over-clustering: a shard that is a whole population | Over-clustering control: cluster_2 is a distinct population (PBMC3k) | ROBUST | **ROBUST** | 23/23 | 100% | control OK |
| 14 | A QC threshold that silently deletes a real population | QC deletion: NK cells (GNLY+ NKG7+) exist (PBMC3k) | not ROBUST; the failure table names max_mito_pct | **NOT SUPPORTED** | 4/23 | 17% | caught |
| 15 | A fake batch effect confounded with a subset of cells | Batch artifact: 'CD4 T subtype' is a distinct population (PBMC3k) | not ROBUST; holds only when batch correction is off | **NOT SUPPORTED** | 4/24 | 17% | caught |
| 16 | Control for the batch effect: the real parent population | Batch control: 'CD4 T' is a distinct population (PBMC3k) | ROBUST | **ROBUST** | 24/24 | 100% | control OK |
