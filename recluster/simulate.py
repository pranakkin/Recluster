"""Synthetic single-cell counts with known structure, for the tests and the synthetic example.

This is deliberately simple (gamma-Poisson counts, marker genes with fold changes, donor
random effects). It is *not* used to produce any benchmark headline; the real-data benchmark
uses Kang et al. 2018. Its job is to give the test-suite claims that must pass and claims
that must fail by construction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


@dataclass
class SimConfig:
    n_donors: int = 6
    conditions: tuple[str, str] = ("ctrl", "stim")
    cells_per_sample: int = 150
    n_genes: int = 2500
    n_markers_per_type: int = 15
    marker_fold: float = 6.0
    cell_types: tuple[str, ...] = ("T", "B", "Mono", "NK")
    base_proportions: tuple[float, ...] = (0.45, 0.20, 0.25, 0.10)
    # planted composition shift: Mono more abundant in stim
    shift_type: str | None = "Mono"
    shift_proportion_in_stim: float = 0.40
    # planted DE: genes up in Mono under stim
    de_type: str | None = "Mono"
    n_de_genes: int = 5
    de_fold: float = 4.0
    donor_sd: float = 0.15
    dispersion: float = 0.3
    n_mito: int = 10
    mito_fraction: float = 0.05
    seed: int = 0


def simulate(cfg: SimConfig | None = None) -> AnnData:
    """Simulate counts. Truth is stored in ``.obs`` and ``.uns['truth']``."""
    cfg = cfg or SimConfig()
    rng = np.random.default_rng(cfg.seed)
    types = list(cfg.cell_types)
    n_types = len(types)
    n_genes = cfg.n_genes

    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    # gene-level baseline means (log-normal), heavy tail like real data
    base_mean = np.exp(rng.normal(-1.0, 1.2, size=n_genes)).astype(np.float64)

    # marker blocks: disjoint sets of genes per type
    marker_idx: dict[str, np.ndarray] = {}
    pool = rng.permutation(n_genes)
    cursor = 0
    for t in types:
        marker_idx[t] = pool[cursor : cursor + cfg.n_markers_per_type]
        cursor += cfg.n_markers_per_type
        for j, g in enumerate(marker_idx[t]):
            gene_names[g] = f"{t}_marker_{j}"

    # DE genes (disjoint from markers)
    de_idx = pool[cursor : cursor + cfg.n_de_genes]
    cursor += cfg.n_de_genes
    for j, g in enumerate(de_idx):
        gene_names[g] = f"DE_UP_{j}"
    # null genes with a name that makes the "should fail" tests readable
    null_idx = pool[cursor : cursor + cfg.n_de_genes]
    cursor += cfg.n_de_genes
    for j, g in enumerate(null_idx):
        gene_names[g] = f"NULL_{j}"
    mito_idx = pool[cursor : cursor + cfg.n_mito]
    for j, g in enumerate(mito_idx):
        gene_names[g] = f"MT-{j}"

    # per-type expression profile
    profiles = np.tile(base_mean, (n_types, 1))
    for k, t in enumerate(types):
        profiles[k, marker_idx[t]] *= cfg.marker_fold
        # make markers reasonably expressed
        profiles[k, marker_idx[t]] = np.maximum(profiles[k, marker_idx[t]], 1.5)
    for k in range(n_types):
        profiles[k, mito_idx] = base_mean.sum() * cfg.mito_fraction / cfg.n_mito
    for k in range(n_types):
        profiles[k, de_idx] = np.maximum(profiles[k, de_idx], 1.0)
        profiles[k, null_idx] = np.maximum(profiles[k, null_idx], 1.0)

    obs_rows: list[dict[str, object]] = []
    blocks: list[sp.csr_matrix] = []
    base_props = np.asarray(cfg.base_proportions, dtype=float)
    base_props = base_props / base_props.sum()

    for d in range(cfg.n_donors):
        donor = f"donor{d}"
        donor_effect = np.exp(rng.normal(0, cfg.donor_sd, size=n_genes))
        for cond in cfg.conditions:
            props = base_props.copy()
            if cfg.shift_type is not None and cond == cfg.conditions[1]:
                k = types.index(cfg.shift_type)
                rest = 1.0 - cfg.shift_proportion_in_stim
                others = np.delete(props, k)
                props[k] = cfg.shift_proportion_in_stim
                props[np.arange(n_types) != k] = others / others.sum() * rest
            # sample-level noise on proportions
            noisy = rng.dirichlet(props * 60)
            counts_per_type = rng.multinomial(cfg.cells_per_sample, noisy)
            for k, t in enumerate(types):
                n = int(counts_per_type[k])
                if n == 0:
                    continue
                mu = profiles[k] * donor_effect
                if cfg.de_type == t and cond == cfg.conditions[1]:
                    mu = mu.copy()
                    mu[de_idx] *= cfg.de_fold
                size_factor = np.exp(rng.normal(0, 0.25, size=n))[:, None]
                lam = mu[None, :] * size_factor
                # gamma-Poisson for over-dispersion
                shape = 1.0 / cfg.dispersion
                lam = rng.gamma(shape, lam / shape)
                x = rng.poisson(lam).astype(np.float32)
                blocks.append(sp.csr_matrix(x))
                for _ in range(n):
                    obs_rows.append(
                        {
                            "cell_type": t,
                            "donor": donor,
                            "condition": cond,
                            "sample": f"{donor}|{cond}",
                        }
                    )

    X = sp.vstack(blocks).tocsr()
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell{i}" for i in range(len(obs))]
    for c in ("cell_type", "donor", "condition", "sample"):
        obs[c] = obs[c].astype("category")
    adata = AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))
    adata.uns["truth"] = {
        "markers": {t: [gene_names[g] for g in marker_idx[t]] for t in types},
        "de_genes": [gene_names[g] for g in de_idx],
        "null_genes": [gene_names[g] for g in null_idx],
        "shift_type": cfg.shift_type,
        "de_type": cfg.de_type,
    }
    return adata
