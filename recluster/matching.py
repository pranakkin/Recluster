"""Find the claimed population inside one grid run.

Signals: marker enrichment (the stratified within-sample AUROC of a per-cell marker score for a
cluster versus the other cells of the same sample, pooled over samples with weight n1 * n0, so
that a sample-level shift of the score or a QC setting that changes the mix of other cell types
in one sample cannot move the other samples' terms), per-marker enrichment (a marker set must be
carried by most of its genes, not one), label overlap (fraction of a cluster carrying the claimed
label), and overlap with a reference cell set (Jaccard over barcodes that passed the run's QC,
with the QC-removed fraction reported separately). A cluster can qualify by markers, by labels,
or by either (``Thresholds.candidate_rule``). Qualifying by labels needs a reference built from
the claim's label column; a reference carried over from the baseline run's match is only
compared with.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import rankdata

from recluster.config import Thresholds


def lognorm_matrix(counts: AnnData, genes: list[str]) -> pd.DataFrame:
    """CP10k log1p expression for a few genes, all cells. Independent of any grid choice."""
    if not genes:
        return pd.DataFrame(index=counts.obs_names)
    x = counts.X
    total = np.asarray(x.sum(axis=1)).ravel()
    idx = counts.var_names.get_indexer(genes)
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise KeyError(f"genes not in counts matrix: {missing}")
    sub = x[:, idx]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.where(total[:, None] > 0, sub / total[:, None] * 1e4, 0.0)
    return pd.DataFrame(np.log1p(norm), index=counts.obs_names, columns=genes)


def _zscore(s: pd.Series) -> pd.Series:
    sd = float(s.std(ddof=0))
    return (s - s.mean()) / sd if sd > 0 else s * 0.0


def marker_score(
    lognorm: pd.DataFrame, markers: list[str], negative: list[str], logic: str = "any"
) -> pd.Series:
    """Per-cell marker score, z-scored over cells.

    ``any``: mean marker expression minus mean negative-marker expression (markers are
    alternative evidence for one cell type; tolerant to dropout).
    ``all``: the minimum of the per-marker z-scores (the population must express every
    marker; the natural reading of "CD14+ MS4A1+").
    """
    if logic == "all":
        s = pd.concat([_zscore(lognorm[g]) for g in markers], axis=1).min(axis=1)
    else:
        s = lognorm[markers].mean(axis=1)
    if negative:
        s = s - lognorm[negative].mean(axis=1)
    return _zscore(s).rename("marker_score")


def auroc_by_group(score: np.ndarray, groups: np.ndarray) -> dict[str, float]:
    """AUROC of ``score`` for each group vs the rest, via the Mann-Whitney U identity."""
    ranks = rankdata(score)
    n = len(score)
    out: dict[str, float] = {}
    for g in np.unique(groups):
        mask = groups == g
        n1 = int(mask.sum())
        n0 = n - n1
        if n1 == 0 or n0 == 0:
            out[str(g)] = float("nan")
            continue
        u = ranks[mask].sum() - n1 * (n1 + 1) / 2.0
        out[str(g)] = float(u / (n1 * n0))
    return out


def auroc_two_sets(score: np.ndarray, positive: np.ndarray) -> float:
    """AUROC of ``score`` for the boolean ``positive`` mask vs the rest."""
    n1 = int(positive.sum())
    n0 = len(score) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(score)
    u = ranks[positive].sum() - n1 * (n1 + 1) / 2.0
    return float(u / (n1 * n0))


def within_sample_auroc(
    score: np.ndarray, groups: np.ndarray, samples: np.ndarray
) -> dict[str, float]:
    """Stratified within-sample AUROC of ``score`` for each cluster vs the rest of its samples.

    Each sample contributes the Mann-Whitney U of the cluster's cells against the other cells of
    the same sample, and the AUROC is the summed U over the summed number of (cluster, other)
    cell pairs: a sample's weight is n1 * n0, so tiny samples contribute little and no single
    sample decides. Cells are only ever ranked against cells of their own sample, so a
    sample-level shift of the score, or a QC setting that changes the mix of other cell types in
    one sample, cannot move the other samples' terms, while the whole-dataset AUROC
    (``auroc_by_group``) moves with both. Taking the worst sample instead would reject a rare
    population whenever one sample's few cells happen to rank low. Every sample holding at least
    one cluster cell and at least one other cell contributes; clusters with no contributing
    sample are absent from the result.
    """
    u_sum: dict[str, float] = {}
    pairs: dict[str, int] = {}
    for s in np.unique(samples):
        in_s = samples == s
        ranks = rankdata(score[in_s])
        g_s = groups[in_s]
        n = int(in_s.sum())
        for c in np.unique(g_s):
            mask = g_s == c
            n1 = int(mask.sum())
            n0 = n - n1
            if n0 == 0:
                continue
            key = str(c)
            u_sum[key] = u_sum.get(key, 0.0) + float(ranks[mask].sum() - n1 * (n1 + 1) / 2.0)
            pairs[key] = pairs.get(key, 0) + n1 * n0
    return {c: u_sum[c] / pairs[c] for c in u_sum}


def per_gene_aurocs(values: pd.DataFrame, groups: np.ndarray) -> pd.DataFrame:
    """AUROC of every column (gene) for every group vs the rest: groups x genes."""
    ranks = np.column_stack([rankdata(values[c].to_numpy()) for c in values.columns])
    n = len(groups)
    rows = {}
    for g in np.unique(groups):
        mask = groups == g
        n1 = int(mask.sum())
        n0 = n - n1
        if n1 == 0 or n0 == 0:
            rows[str(g)] = np.full(values.shape[1], np.nan)
            continue
        u = ranks[mask].sum(axis=0) - n1 * (n1 + 1) / 2.0
        rows[str(g)] = u / (n1 * n0)
    return pd.DataFrame(rows, index=values.columns).T


@dataclass
class Match:
    """Where the claimed population landed in one run.

    ``candidate_by`` says how each candidate cluster qualified ("markers", "labels" or "both").
    ``judged_by`` says which statistic judged each marker candidate: "within-sample" (the
    stratified within-sample AUROC) or "overall" (the whole-dataset AUROC, used for a cluster
    with no contributing sample and when there is no sample column). ``marker_basis`` sums that
    up for the run as "within-sample", "overall" or "mixed"; with no marker candidate it names
    the statistic behind ``best_cluster_auroc``.
    """

    mode: str  # "markers", "labels" or "either"
    clusters: list[str]
    n_cells: int
    found: bool
    best_cluster_auroc: float | None  # the best judged value over all clusters
    marker_basis: str | None  # "within-sample", "overall" or "mixed"
    marker_auroc: float | None
    per_marker_auroc: dict[str, float] | None
    markers_enriched: str | None  # e.g. "3/3"
    jaccard: float | None
    precision: float | None
    recall: float | None
    n_reference: int | None
    reference_qc_removed: float | None  # fraction of reference cells this run's QC removed
    n_samples_present: int | None
    max_sample_share: float | None
    candidate_by: dict[str, str] = field(default_factory=dict)
    judged_by: dict[str, str] = field(default_factory=dict)
    cells: pd.Index = field(default_factory=lambda: pd.Index([]))

    def to_dict(self) -> dict[str, object]:
        d = asdict(self)
        d.pop("cells")
        return d


def _marker_candidates(
    clusters_arr: np.ndarray,
    sizes: pd.Series,
    score: np.ndarray,
    marker_frame: pd.DataFrame | None,
    logic: str,
    thr: Thresholds,
    samples: np.ndarray | None,
) -> tuple[list[str], float | None, dict[str, str], str]:
    """Clusters that qualify by markers.

    Returns the candidates, the best judged AUROC over all clusters, ``judged_by`` for the
    candidates and the run's ``marker_basis`` (see ``Match``). A cluster is judged by its
    stratified within-sample AUROC when any sample contributes one, otherwise by the
    whole-dataset AUROC.
    """
    overall = auroc_by_group(score, clusters_arr)
    within = within_sample_auroc(score, clusters_arr, samples) if samples is not None else {}
    judged = {c: within.get(c, a) for c, a in overall.items()}
    statistic = {c: "within-sample" if c in within else "overall" for c in overall}
    finite = {c: a for c, a in judged.items() if not np.isnan(a)}
    best_cluster = max(finite, key=finite.__getitem__) if finite else None
    best = finite[best_cluster] if best_cluster is not None else None
    per_gene = None
    if marker_frame is not None and marker_frame.shape[1] > 1:
        per_gene = per_gene_aurocs(marker_frame, clusters_arr)
    out = []
    for c, a in judged.items():
        if not (a >= thr.marker_auroc_min and sizes[c] >= thr.min_cluster_size):
            continue
        if per_gene is not None:
            enriched = (per_gene.loc[c] >= thr.marker_per_gene_auroc_min).to_numpy()
            if logic == "all" and not enriched.all():
                continue
            if logic != "all" and not (enriched.mean() > thr.marker_majority_min):
                continue
        out.append(c)
    judged_by = {c: statistic[c] for c in out}
    kinds = set(judged_by.values())
    if len(kinds) > 1:
        basis = "mixed"
    elif kinds:
        basis = kinds.pop()
    else:
        basis = statistic[best_cluster] if best_cluster is not None else "overall"
    return out, best, judged_by, basis


def match_population(
    clusters: pd.Series,
    *,
    score: pd.Series | None,
    reference: pd.Index | None,
    sample: pd.Series | None,
    min_cells: int,
    thr: Thresholds,
    marker_frame: pd.DataFrame | None = None,
    logic: str = "any",
    reference_is_labelled: bool = False,
) -> Match:
    """Locate the population in a run by marker enrichment and/or label overlap.

    A cluster qualifies by markers when its stratified within-sample AUROC (the whole-dataset
    AUROC without a sample column, or for a cluster no sample can judge) reaches
    ``thr.marker_auroc_min`` and most of its markers are individually enriched, and by labels
    when at least ``thr.label_precision_min`` of its cells are in ``reference``. Qualifying by
    labels needs ``reference_is_labelled``: a reference built from the claim's label column
    records which cells carry the label, whereas one carried over from an earlier run's match is
    the thing under test, so it is used only for the Jaccard/precision/recall overlap. Without
    this rule a merged low-resolution cluster holding most of the baseline population would
    qualify with no marker evidence at all.
    """
    run_cells = clusters.index
    labels = clusters.to_numpy().astype(str)
    sizes = pd.Series(labels).value_counts()
    best_auroc: float | None = None
    basis: str | None = None
    judged_by: dict[str, str] = {}
    by_markers: list[str] = []
    by_labels: list[str] = []

    labelled = reference is not None and reference_is_labelled
    if score is None and not labelled:
        raise ValueError("matching needs a marker score and/or a labelled reference set")
    use_markers = score is not None and (thr.candidate_rule != "labels" or not labelled)
    use_labels = labelled and (thr.candidate_rule != "markers" or score is None)

    if use_markers and score is not None:
        s = score.loc[run_cells].to_numpy()
        frame = marker_frame.loc[run_cells] if marker_frame is not None else None
        samples = sample.loc[run_cells].astype(str).to_numpy() if sample is not None else None
        by_markers, best_auroc, judged_by, basis = _marker_candidates(
            labels, sizes, s, frame, logic, thr, samples
        )
    if use_labels and reference is not None:
        in_ref = pd.Series(run_cells.isin(reference), index=run_cells)
        precision_by_cluster = in_ref.groupby(labels).mean()
        by_labels = [
            str(c)
            for c, p in precision_by_cluster.items()
            if p >= thr.label_precision_min and sizes[str(c)] >= thr.min_cluster_size
        ]
    candidate_by = {c: "markers" for c in by_markers}
    for c in by_labels:
        candidate_by[c] = "both" if c in candidate_by else "labels"
    candidates = list(candidate_by)
    mode = "either" if (use_markers and use_labels) else ("markers" if use_markers else "labels")

    in_m = np.isin(labels, candidates)
    cells = run_cells[in_m]
    n_cells = int(len(cells))

    marker_auroc = None
    per_marker: dict[str, float] | None = None
    enriched_str = None
    if score is not None and 0 < n_cells < len(run_cells):
        marker_auroc = auroc_two_sets(score.loc[run_cells].to_numpy(), in_m)
        if marker_frame is not None:
            group = np.where(in_m, "pop", "rest")
            pg = per_gene_aurocs(marker_frame.loc[run_cells], group)
            if "pop" in pg.index:
                per_marker = {g: float(v) for g, v in pg.loc["pop"].items()}
                n_enr = int(sum(v >= thr.marker_per_gene_auroc_min for v in per_marker.values()))
                enriched_str = f"{n_enr}/{len(per_marker)}"

    jaccard = precision = recall = None
    n_reference = None
    ref_removed = None
    if reference is not None:
        n_reference = int(len(reference))
        ref_kept = reference.intersection(run_cells) if thr.jaccard_on_kept_cells else reference
        ref_removed = 1.0 - len(ref_kept) / n_reference if n_reference else None
        inter = int(cells.isin(ref_kept).sum())
        union = n_cells + len(ref_kept) - inter
        jaccard = inter / union if union else 0.0
        precision = inter / n_cells if n_cells else 0.0
        recall = inter / len(ref_kept) if len(ref_kept) else 0.0

    n_samples_present = max_share = None
    if sample is not None and n_cells > 0:
        per_sample = sample.loc[cells].value_counts()
        n_samples_present = int((per_sample >= thr.min_cells_per_sample).sum())
        max_share = float(per_sample.max() / n_cells)
    elif sample is not None:
        n_samples_present, max_share = 0, None

    return Match(
        mode=mode,
        clusters=sorted(candidates, key=_cluster_sort_key),
        n_cells=n_cells,
        found=n_cells >= min_cells,
        best_cluster_auroc=best_auroc,
        marker_basis=basis,
        marker_auroc=marker_auroc,
        per_marker_auroc=per_marker,
        markers_enriched=enriched_str,
        jaccard=jaccard,
        precision=precision,
        recall=recall,
        n_reference=n_reference,
        reference_qc_removed=ref_removed,
        n_samples_present=n_samples_present,
        max_sample_share=max_share,
        candidate_by={c: candidate_by[c] for c in sorted(candidate_by, key=_cluster_sort_key)},
        judged_by={c: judged_by[c] for c in sorted(judged_by, key=_cluster_sort_key)},
        cells=cells,
    )


def _cluster_sort_key(c: str) -> tuple[int, str]:
    return (int(c), "") if c.isdigit() else (10**9, c)
