"""population_exists: is the population recovered as a distinct cluster, across samples?"""

from __future__ import annotations

from recluster.config import Thresholds
from recluster.matching import Match


def evaluate_population(match: Match, thr: Thresholds) -> dict[str, object]:
    """Per-run outcome for a population_exists claim.

    recovered   : enough cells, marker AUROC and Jaccard above their thresholds
    replicate_ok: present in at least ``min_samples_present`` samples (when samples exist)
    holds       : recovered and replicate_ok
    """
    auroc_ok = match.marker_auroc is None or match.marker_auroc >= thr.marker_auroc_min
    jaccard_ok = match.jaccard is None or match.jaccard >= thr.jaccard_min
    recovered = bool(match.found and auroc_ok and jaccard_ok)
    replicate_ok = (
        True
        if match.n_samples_present is None
        else match.n_samples_present >= thr.min_samples_present
    )
    reasons: list[str] = []
    if not match.found:
        reasons.append(f"only {match.n_cells} matched cells")
    if not auroc_ok:
        reasons.append(f"marker AUROC {match.marker_auroc:.2f} < {thr.marker_auroc_min}")
    if not jaccard_ok and match.jaccard is not None:
        reasons.append(f"Jaccard {match.jaccard:.2f} < {thr.jaccard_min}")
    if not replicate_ok:
        reasons.append(
            f"present in {match.n_samples_present} sample(s) < {thr.min_samples_present}"
        )
    return {
        "recovered": recovered,
        "replicate_ok": bool(replicate_ok),
        "holds": bool(recovered and replicate_ok),
        "holds_cell_level": recovered,  # existence without the replicate check
        "reasons": reasons,
        "key_stats": {"jaccard": match.jaccard, "marker_auroc": match.marker_auroc},
    }
