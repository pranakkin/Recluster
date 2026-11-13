Scored claims: 40

### pooled

| caller | sensitivity (95 % CI) | specificity (95 % CI) |
|---|---|---|
| Recluster ROBUST | 20/20 (100 %, 84-100) | 17/20 (85 %, 64-95) |
| Recluster ROBUST or FRAGILE | 20/20 (100 %, 84-100) | 17/20 (85 %, 64-95) |
| naive cell-level test (baseline run) | 20/20 (100 %, 84-100) | 0/20 (0 %, 0-16) |
| one pseudobulk analysis (baseline run) | 20/20 (100 %, 84-100) | 17/20 (85 %, 64-95) |

### mouse_lps4

| caller | sensitivity (95 % CI) | specificity (95 % CI) |
|---|---|---|
| Recluster ROBUST | 10/10 (100 %, 72-100) | 10/10 (100 %, 72-100) |
| Recluster ROBUST or FRAGILE | 10/10 (100 %, 72-100) | 10/10 (100 %, 72-100) |
| naive cell-level test (baseline run) | 10/10 (100 %, 72-100) | 0/10 (0 %, 0-28) |
| one pseudobulk analysis (baseline run) | 10/10 (100 %, 72-100) | 10/10 (100 %, 72-100) |

### rat_lps4

| caller | sensitivity (95 % CI) | specificity (95 % CI) |
|---|---|---|
| Recluster ROBUST | 10/10 (100 %, 72-100) | 7/10 (70 %, 40-89) |
| Recluster ROBUST or FRAGILE | 10/10 (100 %, 72-100) | 7/10 (70 %, 40-89) |
| naive cell-level test (baseline run) | 10/10 (100 %, 72-100) | 0/10 (0 %, 0-28) |
| one pseudobulk analysis (baseline run) | 10/10 (100 %, 72-100) | 7/10 (70 %, 40-89) |

### Misses and false calls

- Recluster ROBUST: misses none; false calls ['rat_lps4/Alkbh7', 'rat_lps4/Ezh2', 'rat_lps4/Son']
- Recluster ROBUST or FRAGILE: misses none; false calls ['rat_lps4/Alkbh7', 'rat_lps4/Ezh2', 'rat_lps4/Son']
- naive cell-level test (baseline run): misses none; false calls ['mouse_lps4/2210016L21Rik', 'mouse_lps4/Anapc5', 'mouse_lps4/Fam193a', 'mouse_lps4/Itgb7', 'mouse_lps4/Nmt1', 'mouse_lps4/Rab14', 'mouse_lps4/Selenok', 'mouse_lps4/Snrpa', 'mouse_lps4/Tbl1x', 'mouse_lps4/Uimc1', 'rat_lps4/Aarsd1', 'rat_lps4/Alkbh7', 'rat_lps4/Ezh2', 'rat_lps4/Ing3', 'rat_lps4/Nipsnap3b', 'rat_lps4/Prpf6', 'rat_lps4/Serinc2', 'rat_lps4/Son', 'rat_lps4/Tbl3', 'rat_lps4/Ubap2']
- one pseudobulk analysis (baseline run): misses none; false calls ['rat_lps4/Alkbh7', 'rat_lps4/Ezh2', 'rat_lps4/Son']

### Every claim

| claim | truth | bulk log2FC | verdict | runs | naive cell | one pseudobulk |
|---|---|---|---|---|---|---|
| mouse_lps4/Acsl1 | change | 4.23 | ROBUST | 23/23 | call | call |
| mouse_lps4/Ampd3 | change | 3.01 | ROBUST | 23/23 | call | call |
| mouse_lps4/Azi2 | change | 1.74 | ROBUST | 23/23 | call | call |
| mouse_lps4/Cry1 | change | 1.25 | ROBUST | 23/23 | call | call |
| mouse_lps4/Flvcr1 | change | -1.95 | ROBUST | 23/23 | call | call |
| mouse_lps4/Gprc5b | change | -2.23 | ROBUST | 23/23 | call | call |
| mouse_lps4/Mtmr14 | change | 2.31 | ROBUST | 23/23 | call | call |
| mouse_lps4/Ppfia3 | change | 4.69 | ROBUST | 23/23 | call | call |
| mouse_lps4/Spred1 | change | 3.40 | ROBUST | 23/23 | call | call |
| mouse_lps4/Tap1 | change | 3.36 | ROBUST | 23/23 | call | call |
| mouse_lps4/2210016L21Rik | no_change | -0.09 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Anapc5 | no_change | -0.16 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Fam193a | no_change | -0.16 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Itgb7 | no_change | -0.28 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Nmt1 | no_change | 0.06 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Rab14 | no_change | 0.10 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Selenok | no_change | 0.16 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Snrpa | no_change | -0.21 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Tbl1x | no_change | 0.05 | NOT SUPPORTED | 0/23 | call | - |
| mouse_lps4/Uimc1 | no_change | 0.01 | NOT SUPPORTED | 0/23 | call | - |
| rat_lps4/Ada | change | 1.53 | ROBUST | 21/21 | call | call |
| rat_lps4/Ankrd33b | change | 1.43 | ROBUST | 21/21 | call | call |
| rat_lps4/Baz1b | change | -1.01 | ROBUST | 21/21 | call | call |
| rat_lps4/Cmpk2 | change | 6.60 | ROBUST | 21/21 | call | call |
| rat_lps4/Etv1 | change | -1.89 | ROBUST | 21/21 | call | call |
| rat_lps4/Fnip2 | change | -1.76 | ROBUST | 21/21 | call | call |
| rat_lps4/Mtmr10 | change | -1.03 | ROBUST | 21/21 | call | call |
| rat_lps4/Polr1e | change | 2.53 | ROBUST | 21/21 | call | call |
| rat_lps4/Sp110 | change | 1.06 | ROBUST | 21/21 | call | call |
| rat_lps4/Tapbp | change | 1.37 | ROBUST | 21/21 | call | call |
| rat_lps4/Aarsd1 | no_change | -0.08 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Alkbh7 | no_change | 0.10 | ROBUST | 21/21 | call | call |
| rat_lps4/Ezh2 | no_change | -0.32 | ROBUST | 21/21 | call | call |
| rat_lps4/Ing3 | no_change | 0.13 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Nipsnap3b | no_change | -0.04 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Prpf6 | no_change | 0.08 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Serinc2 | no_change | -0.04 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Son | no_change | -0.03 | ROBUST | 21/21 | call | call |
| rat_lps4/Tbl3 | no_change | 0.09 | NOT SUPPORTED | 0/21 | call | - |
| rat_lps4/Ubap2 | no_change | 0.05 | NOT SUPPORTED | 0/21 | call | - |
