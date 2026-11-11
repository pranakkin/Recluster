# Benchmark

Two parts, one command. Results exactly as measured are in [RESULTS.md](RESULTS.md); the
machine-readable version is `results.json`.

```bash
uv run recluster bench
```

That downloads what it needs into `data/` (gitignored), builds every case deterministically
(seed 0), runs each claim through the normal engine with the default grid and thresholds, and
writes `bench/results.json`, `bench/results_table.md` and one report folder per claim under
`bench/out/`. Use `--part kang` or `--part planted` for one half and `--jobs N` to change how
many claims run in parallel (default 4). The whole thing takes about 40 minutes on a laptop.

## Part 1: real data, real ground truth (Kang et al. 2018)

PBMCs from 8 lupus patients, each split into an untreated and a 6 h IFN-beta-stimulated sample
and demultiplexed by genotype. Source: the processed object distributed with pertpy
(`https://exampledata.scverse.org/pertpy/kang_2018.h5ad`, SHA-256 verified), which carries raw
counts, `label` (ctrl/stim), `replicate` (donor) and `cell_type`.

- **True claims**, expected ROBUST: interferon-stimulated genes up in monocytes and T cells;
  monocytes and NK cells present across donors.
- **False claims built from the real data**, expected NOT SUPPORTED: control cells only with
  donors split into two fake groups, and all cells with condition labels shuffled at the
  donor level. The cell-level column shows whether naive statistics would have called them.
- **Doublets**: the processed object has singlets only, so Scrublet is compared with demuxlet's
  calls on the two original 10x lanes from GEO (GSE96583, batch 2: 29,065 barcodes).

## Part 2: planted errors (PBMC3k)

Each error is built by a fixed function of `scanpy.datasets.pbmc3k()` and seed 0, and its
"caught" rule is written in `recluster/bench/planted.py` before the run:

1. synthetic doublets (CD14 monocyte + B cell) appended as a fake population;
2. over-clustering at Leiden resolution 3.0, claiming one shard of CD4 T cells (with a control);
3. the mitochondrial fraction of NK cells raised into 12 to 16 %, so the baseline QC deletes them;
4. a gene-level batch effect applied to a random 40 % of cells, claiming a batch-pure cluster
   as a subtype (with a control).

Thresholds and grid are the package defaults. Any change made after seeing results is logged in
`RESULTS.md` with the reason.
