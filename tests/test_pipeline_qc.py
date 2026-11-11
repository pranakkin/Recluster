"""QC-cap bookkeeping: capped values print cleanly, grid points that share a cached stage warn
about their own requested value, effective params fold the cap in, and the Scrublet fallback
quotes its rate as a percentage."""

from __future__ import annotations

import re

import numpy as np
import pytest
from anndata import AnnData

from recluster.config import GridConfig, PipelineParams
from recluster.data import detect_mito_prefix, extract_counts
from recluster.pipeline import Pipeline

CAPPED_TO = re.compile(r"capped to (\S+)")
ONE_DECIMAL = re.compile(r"\d+(\.\d)?")


@pytest.fixture(scope="module")
def pipe(sim: AnnData) -> Pipeline:
    counts, _ = extract_counts(sim)
    return Pipeline(counts, GridConfig(), mito_prefix=detect_mito_prefix(counts.var_names))


def _capping_mito(pipe: Pipeline, keep_fraction: float) -> float:
    """A max_mito_pct that keeps about ``keep_fraction`` of the cells, so the cap pulls it back."""
    pct = pipe.counts.obs["pct_counts_mt"].to_numpy()
    value = round(float(np.quantile(pct, keep_fraction)), 1)
    assert (pct > value).mean() > pipe.cfg.qc_max_cell_loss
    return value


def test_capped_max_mito_has_one_decimal_and_prints_cleanly(pipe: Pipeline) -> None:
    harsh = _capping_mito(pipe, 0.5)
    min_genes, max_mito, capped = pipe.effective_qc(PipelineParams(max_mito_pct=harsh))
    assert capped == {"max_mito_pct": max_mito} and max_mito > harsh
    assert type(max_mito) is float and max_mito == round(max_mito, 1)
    assert isinstance(min_genes, int)
    run = pipe.run(PipelineParams(max_mito_pct=harsh, batch_correction=False))
    assert run.effective["max_mito_pct"] == max_mito and run.effective["capped"] == capped
    # the label fragment the engine builds from effective["capped"], plus the warning line
    texts = [f"max_mito_pct capped to {capped['max_mito_pct']}"]
    texts += [w for w in run.warnings if "capped to" in w]
    assert len(texts) == 2
    numbers = CAPPED_TO.findall(" ".join(texts))
    assert len(numbers) == 2
    for number in numbers:
        assert ONE_DECIMAL.fullmatch(number), texts


def test_runs_sharing_a_capped_stage_warn_about_their_own_value(pipe: Pipeline) -> None:
    a, b = _capping_mito(pipe, 0.5), _capping_mito(pipe, 0.4)
    assert a != b
    pa = PipelineParams(max_mito_pct=a, batch_correction=False)
    pb = PipelineParams(max_mito_pct=b, batch_correction=False)
    ea, eb = pipe.effective_params(pa), pipe.effective_params(pb)
    assert ea == eb and ea != pa
    assert ea.max_mito_pct == pipe.effective_qc(pa)[1] and ea.min_genes == pa.min_genes
    assert pipe.qc_stage(pa).adata is pipe.qc_stage(pb).adata  # one cached stage for both
    loss = f"{pipe.cfg.qc_max_cell_loss:.0%}"
    for params in (pa, pb, pa):
        expected = (
            f"max_mito_pct={params.max_mito_pct} would remove more than {loss} of cells; "
            f"capped to {ea.max_mito_pct}"
        )
        for warns in (pipe.qc_stage(params).warnings, pipe.run(params).warnings):
            assert [w for w in warns if "capped to" in w] == [expected]


def test_effective_params_is_identity_without_capping(pipe: Pipeline, sim: AnnData) -> None:
    p = PipelineParams(min_genes=50, max_mito_pct=20.0, resolution=0.5)
    assert pipe.effective_qc(p)[2] == {}
    assert pipe.effective_params(p) == p
    # without mitochondrial genes the mito axis is never capped
    counts, _ = extract_counts(sim)
    no_mito = Pipeline(counts, GridConfig(), mito_prefix=None)
    harsh = PipelineParams(max_mito_pct=0.5)
    assert no_mito.effective_params(harsh) == harsh


def test_fallback_warning_quotes_rate_as_percentage(
    monkeypatch: pytest.MonkeyPatch, sim: AnnData
) -> None:
    import scanpy as sc

    counts, _ = extract_counts(sim)
    pipe = Pipeline(counts, GridConfig(), mito_prefix="MT-")

    def fake_scrublet(adata: AnnData, **kw: object) -> None:  # calls nothing: fallback is used
        rng = np.random.default_rng(0)
        adata.obs["doublet_score"] = rng.uniform(0, 0.3, adata.n_obs)
        adata.obs["predicted_doublet"] = False

    monkeypatch.setattr(sc.pp, "scrublet", fake_scrublet)
    stage = pipe.qc_stage(PipelineParams())
    assert stage.effective["scrublet_fallback_subsets"] == 1
    line = next(w for w in stage.warnings if "expected rate" in w)
    assert "expected rate 0.8% per 1,000 cells" in line and "0.008" not in line
