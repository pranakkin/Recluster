# Benchmark

Two parts, one command. Results are in [RESULTS.md](RESULTS.md); the machine-readable version
is `results.json`.

```bash
uv run recluster bench
```

That downloads what it needs into `data/` (gitignored), builds every case deterministically
(fixed seeds), runs each claim through the normal engine with the default grid and thresholds, and
writes `bench/results.json`, `bench/results_table.md` and one report folder per claim under
`bench/out/`. Use `--part kang` or `--part planted` for one half and `--jobs N` to change how
many claims run in parallel (default 4). The whole benchmark took 38 minutes with 4 workers on
an Apple M5 Pro laptop. The run rewrites the committed `bench/results.json` and
`bench/results_table.md`; with `--part`, the rewritten files hold only that part's cases.
`git diff bench/` shows any difference and
`git checkout bench/results.json bench/results_table.md` restores them.

## Part 1: real data (Kang et al. 2018)

PBMCs from 8 lupus patients, each split into an untreated and a 6 h IFN-beta-stimulated sample
and demultiplexed by genotype. Source: the processed object distributed with pertpy
(`https://exampledata.scverse.org/pertpy/kang_2018.h5ad`, SHA-256 verified), which carries raw
counts, `label` (ctrl/stim), `replicate` (donor) and `cell_type`.

- **True claims**, expected ROBUST: interferon-stimulated genes up in monocytes and T cells;
  monocytes and NK cells present across donors.
- **False claims built from the real data**, expected NOT SUPPORTED: control cells only with
  donors split into two fake groups, and all cells with condition labels shuffled at the
  donor level. The cell-level column shows whether naive statistics would have called them.
- **A planted time-of-day effect**: the 8 donors split into a morning group and an afternoon
  group (one collection time per donor), with a documented clock-gene shift and a rhythmic gene
  (SQSTM1) planted into every cell as a function of the hour, then claimed as a condition
  effect. Run once with the collection times hidden (the clock-gene screen must flag it) and
  once with them given (the adjustment must remove it); the true IFN-β claims above report their
  time-check outcome so that false alarms are visible.
- **Doublets**: the processed object has singlets only, so Scrublet is compared with demuxlet's
  calls on the two original 10x lanes from GEO (GSE96583, batch 2: 29,065 barcodes); the first
  run measured Scrublet's automatic threshold alone, the current one includes the rate-based
  fallback.

## Part 2: planted errors (PBMC3k)

Each error is built by a fixed function of `scanpy.datasets.pbmc3k()` and a fixed seed, and its
"caught" rule is in `recluster/bench/planted.py`:

1. synthetic doublets (CD14 monocyte + B cell) appended as a fake population;
2. over-clustering at Leiden resolution 3.0, claiming one shard of CD4 T cells (with a control);
3. the mitochondrial fraction of NK cells raised into 12 to 16 %, so the baseline QC deletes them;
4. a gene-level batch effect applied to a random 40 % of cells, claiming a batch-pure cluster
   as a subtype (with a control);
5. two held-out doublet blocks (NK + CD14 monocyte, CD8 T + B) appended as fake populations.
   They and their rules were written into `PLAN.md` (commit 6201b82) before the
   doublet-population step was built, and were run once.

The rules for errors 1 to 4 were committed together with the first results (6c2a173); the rule
for error 1 was tightened before the second run (`RESULTS.md`, note C).

Thresholds and grid are the package defaults. Any change made after seeing results is logged in
`RESULTS.md` with the reason, and the first run's table is kept next to the current one
(`results_v0.1_table.md`, `results_v0.1.json`).
