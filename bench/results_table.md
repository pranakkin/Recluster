| # | Case | Claim | Expected | Verdict | Held | Cell-level | Time check | Outcome |
|---|---|---|---|---|---|---|---|---|
| 1 | ISG15 up in CD14+ monocytes after IFN-beta | ISG15 up in CD14+ monocytes, stim vs ctrl (Kang 2018) | ROBUST | **ROBUST** | 21/21 | 100% | none | correct |
| 2 | Four interferon-stimulated genes up in CD14+ monocytes | ISG15, IFI6, IFIT3 and MX1 all up in CD14+ monocytes, stim vs ctrl (Kang 2018) | ROBUST | **ROBUST** | 21/21 | 100% | none | correct |
| 3 | ISG15 up in CD4 T cells after IFN-beta | ISG15 up in CD4 T cells, stim vs ctrl (Kang 2018) | ROBUST | **ROBUST** | 21/21 | 100% | none | correct |
| 4 | CD14+ monocytes exist across donors | CD14+ monocytes exist (Kang 2018) | ROBUST | **ROBUST** | 21/21 | – | – | correct |
| 5 | NK cells exist across donors | NK cells exist (Kang 2018) | ROBUST | **ROBUST** | 20/21 | – | – | correct |
| 6 | CD14+ monocyte proportion differs stim vs ctrl | CD14+ monocytes differ in abundance, stim vs ctrl (Kang 2018) | NOT SUPPORTED (no known shift; reported as measured) | **FRAGILE** | 11/21 | 71% | none | FALSE ALARM |
| 7 | False: ISG15 differs between fake donor groups (control cells only) | FALSE: ISG15 differs in CD14+ monocytes between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 1/21 | 100% ⚑ | none | correct |
| 8 | False: LYZ differs between fake donor groups (control cells only) | FALSE: LYZ differs in CD14+ monocytes between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 100% ⚑ | none | correct |
| 9 | False: CD14+ monocyte proportion differs between fake donor groups | FALSE: CD14+ monocyte proportion differs between fake groups A and B (Kang ctrl only) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 24% | none | correct |
| 10 | False: ISG15 up in fake group A (condition shuffled at donor level) | FALSE: ISG15 up in CD14+ monocytes in fake group A vs B (Kang, donor-level shuffle) | NOT SUPPORTED | **NOT SUPPORTED** | 0/21 | 100% ⚑ | none | correct |
| 11 | Planted time-of-day effect claimed as a condition effect (times unknown) | PLANTED TIME, times unknown: SQSTM1 up in CD14+ monocytes, afternoon group B vs morning group A (Kang) | label before the time check ROBUST; time check POSSIBLE TIME CONFOUND; final verdict not ROBUST | **FRAGILE** | 21/21 | 100% | POSSIBLE | caught |
| 12 | Same planted effect with collection times given | PLANTED TIME, times given: SQSTM1 up in CD14+ monocytes, afternoon group B vs morning group A (Kang) | claim held before adjustment; effect removed by the time adjustment; verdict NOT SUPPORTED | **NOT SUPPORTED** | 21/21 | 100% | POSSIBLE | caught |
| 13 | Injected synthetic doublets forming a fake population | Planted doublets: a CD14+ LYZ+ MS4A1+ CD79A+ population exists (PBMC3k) | not ROBUST; holds only when doublet removal is off | **NOT SUPPORTED** | 3/23 | – | – | MISSED |
| 14 | Over-clustering that creates a spurious population | Over-clustering: cluster_14 is a distinct population (PBMC3k) | not ROBUST | **NOT SUPPORTED** | 0/23 | – | – | caught |
| 15 | Control for over-clustering: a shard that is a whole population | Over-clustering control: cluster_2 is a distinct population (PBMC3k) | ROBUST | **ROBUST** | 23/23 | – | – | control OK |
| 16 | A QC threshold that silently deletes a real population | QC deletion: NK cells (GNLY+ NKG7+) exist (PBMC3k) | not ROBUST; the failure table names max_mito_pct | **NOT SUPPORTED** | 4/23 | – | – | caught |
| 17 | A fake batch effect confounded with a subset of cells | Batch artifact: 'CD4 T subtype' is a distinct population (PBMC3k) | not ROBUST; holds only when batch correction is off | **NOT SUPPORTED** | 4/24 | – | – | caught |
| 18 | Control for the batch effect: the real parent population | Batch control: 'CD4 T' is a distinct population (PBMC3k) | ROBUST | **ROBUST** | 24/24 | – | – | control OK |
