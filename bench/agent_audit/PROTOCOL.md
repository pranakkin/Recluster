# Agent audit: protocol

Written before the agent was run, so that the conversion of its claims cannot be fitted to
what it says.

## The run

- Agent: Claude (`claude-opus-5-5`) through the Anthropic Messages API, called by
  `recluster/bench/agent_audit.py`. One tool: a persistent Python 3.11 session in
  `run/workspace/` with the project's installed packages (scanpy, anndata, pydeseq2 and the
  rest). Importing `recluster` is blocked so the agent cannot call the tool that later audits
  it. API keys are removed from the session's environment.
- No system prompt. The agent sees only the user prompt in `run/prompt.md` and the tool
  description in `run/setup.json`. Neither mentions replicates, pseudobulk, donors as the unit
  of analysis, doublets or time of day; the prompt names the donor column only because a user
  would describe their metadata.
- Input: Kang et al. 2018 (GEO GSE96583) as packaged by pertpy, reduced to raw counts plus
  `label` (ctrl/stim) and `replicate` (donor). The authors' cell-type labels, clusters and
  embeddings are removed so the agent does its own annotation.
- Fixed: model id, tool, prompt, `max_tokens` 32000, at most 60 API calls, 20-minute limit per
  code cell. Not fixable: the API has no seed, and this model does not accept `temperature`,
  so a second run would differ. The run is one sample of what the agent does; it is kept
  whole, and its code can be re-executed.
- The run happens once. It is not repeated to obtain a different outcome. If it fails for an
  infrastructure reason (network, API error) before the agent's final answer, that attempt is
  kept in `run_failed_<n>/` and reported, and the run is started again.

## Turning the agent's claims into specs

1. The claims are the agent's final answer (`run/final_answer.md`). Each statement in it that
   asserts something about this dataset is listed in RESULTS.md, quoted exactly.
2. A statement becomes a claim spec when it maps onto a Recluster claim type:
   - "gene G goes up/down in cell type T after stimulation" -> `differential_expression`, one
     spec per gene and cell type, direction as stated (`either` if none is stated);
   - "cell type T is more/less abundant after stimulation" -> `composition_shift`, direction
     as stated (`different` if none is stated);
   - "cell type T exists / is a distinct population" -> `population_exists`.
   Every such statement is tested; none is skipped for being likely to pass or fail. If the
   agent states something that does not map (a pathway, a mechanism, a figure), it is listed as
   not testable, with the reason.
3. The population is the agent's own label: `columns.labels` is the agent's annotation column
   in `annotated.h5ad` and `target.label_value` is the label it used. No markers are added, so
   matching is judged against the agent's labels, as `recluster audit` is meant to be used.
4. Condition `label`, sample `[replicate, label]`, group A `stim`, group B `ctrl`, default grid
   and thresholds. `--jobs` only changes speed.
5. Recluster is run with `recluster audit` on the agent's `annotated.h5ad`. If that file does
   not contain raw counts (`layers['counts']` or `.raw`), the audit input is the original counts
   for the cells the agent kept, with the agent's `.obs` attached by barcode; nothing in the
   agent's annotation is changed, and RESULTS.md says this was needed.
6. The agent's files are never edited. Planted errors, if any are added later, go in a separate,
   labelled section of RESULTS.md and are never mixed with the agent's own claims.

## Amendments (written after reading the agent's answer, before the audits ran)

The agent's answer was read before these were written; no Recluster result on the agent's claims
had been seen except one timing check (amendment 1).

1. **Grouping genes.** Rule 2 said one spec per gene and cell type. The agent names 17 genes as
   going up in all 9 types it tested, so that rule meant about 190 specs. One spec takes 213 s on
   the agent's 24,562 cells with 14 workers (measured on "ISG15 up in CD14 Mono", which came out
   ROBUST, 21/21; that check is one of the agent's claims and is re-run as part of the full set),
   so it would take about 11 hours. Instead, the genes one statement names for one cell type go
   into one spec with `require: all` (the run counts as holding only if every named gene does),
   giving 39 specs. Recluster still records each gene in each run, and RESULTS.md reports any gene
   that fails. One difference from single-gene specs: the pseudobulk p-values are BH-adjusted over
   the genes in the spec, which is stricter than a single-gene spec (where nothing is adjusted).
2. **Groups of cell types.** Where the agent attributes genes to a group ("Monocytes and cDC
   lose ...", "Monocytes and DCs respond most ... Examples are ..."), each gene is tested in each
   member of the group where the agent's own table (`pseudobulk_DE_stim_vs_ctrl.csv`) has it at
   padj < 0.05 in the stated direction, i.e. where the agent's analysis says it. Genes the agent
   pins to one type ("IL27, IDO1 and CH25H in CD16 monocytes", "... in NK") are tested in that
   type only, whatever its table says.
3. **"No change" claims.** The agent says no cell type changes in abundance ("Nothing passes FDR
   < 0.05"). Recluster cannot confirm a null; it tests whether a difference holds. Each of the 11
   types gets a `composition_shift` spec with direction `different`. The agent's statement is
   contradicted for a type only if Recluster finds the shift (ROBUST); NOT SUPPORTED agrees with it,
   and FRAGILE is reported as such.
4. **Named examples only.** The agent says 70 genes rise in all 9 types but names 17. Only the
   named genes are tested: the claim stated in the answer is about those.
5. The specs are produced by `make_specs.py`, which checks that every quote it maps appears
   verbatim in `run/final_answer.md`; `run_audits.py` runs `recluster audit` on each.

### Amendments 6-8 (after an independent cross-check of the mapping, still before any audit)

Three agents extracted the claims from `final_answer.md` without seeing the specs, and a fourth
reviewed `make_specs.py` adversarially. The extractions matched the specs except for one missed
statement; the review found four problems. Changes:

6. **Literal readings are tested too.** Amendment 2 dropped the members of a group for which the
   agent's own table has no support (for example MMP9 and CD93 in CD16 Mono and cDC, where the
   agent never tested them). That favours the agent. These combinations are now tested in specs
   marked *secondary*, reported separately from the primary ones. The missed statement
   ("CXCL10, CCL8 and CXCL11 also appear strongly "up" in T, B and NK cells") is added the same
   way: primary where the agent's table supports it, secondary for the rest.
7. **The abundance null is judged with the agent's own multiplicity.** The agent corrected over
   the 11 types (BH); a Recluster composition spec is tested alone at alpha 0.05, so a shift with
   nominal p < 0.05 would have counted against the agent for that difference alone. Recluster's
   verdicts on the 11 `different` specs are reported as they are, and in addition, within each
   grid run (the 21 runs are the same settings for every claim), Recluster's 11 decisive p-values
   are BH-adjusted together; the agent's "nothing passes FDR < 0.05" is contradicted under the
   grid only if some type passes in at least 80 % of runs (the ROBUST threshold). This
   adjustment is our analysis of Recluster's output, not a Recluster verdict, and RESULTS.md says
   so. A NOT SUPPORTED result is consistent with the agent's null but is not evidence for it
   (8 donor pairs; the agent says the same).
8. **Grouping and `require: all`.** Grouping can turn one weak gene into a failure of the whole
   statement. For every DE spec, each gene's own survival (runs in which that gene holds at the
   replicate level) is reported next to the statement's verdict, from the same runs. The NK
   statement is split into one spec per gene, since two of its genes (AMICA1, NELL2) are not in
   the agent's NK table at all; they come from its CD4 and CD8 T tables, and the sentence still
   says "in NK", so NK is what is tested. Labels the agent hedged (CD4 T, Megakaryocyte/platelet)
   carry the hedge next to the verdict. Statements that map onto no claim type are listed in
   RESULTS.md with the reason.

### Amendment 9 (12 November 2026, after the audits, before the re-run it describes)

9. **Population specs get the sample column.** Checking the results showed that the 11
   population specs (`P-*`) named only the label column. Without a sample column, Recluster runs
   Scrublet on all cells as one batch rather than per donor x condition sample, and at baseline
   it removed 4,826 of 24,562 cells instead of the 349 removed in the comparison specs; up to
   47.7 % of a population (pDC) was lost before matching. That was an error in `make_specs.py`,
   not a choice: the comparison specs already used `sample: [replicate, label]`, and a user
   auditing these labels would give it. The 11 population specs are re-run as `P9-*` with the
   same sample column and nothing else changed. Both sets are reported: the `P-*` verdicts as
   they ran, and the `P9-*` verdicts, which replace them as the primary population results. This
   amendment was written after seeing the `P-*` results, so a reader should weigh any verdict
   that changes with that in mind.

   Erratum to amendment 6: "where the agent never tested them" should read "for which the
   agent's table has no row". The agent's code drops genes both before testing (a count filter)
   and after (rows without an adjusted p), and which applied was not checked.
