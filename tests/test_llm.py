"""The optional LLM layer, tested with a fake client. No network, no key.

What matters: the model only sees column names/values and matching gene symbols; a proposal
that names things not in the data is rejected and retried; summaries cannot introduce numbers
the code did not compute; and everything degrades to a clear message without a key.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from recluster.cli import app
from recluster.llm.common import load_dotenv, numbers_in
from recluster.llm.parse import ProposedClaim, describe_adata, parse_claim
from recluster.llm.summarize import (
    SummaryRejected,
    check_numbers,
    digest_from_provenance,
    summarize_provenance,
)
from recluster.spec import SpecError

from .conftest import mono_markers


class FakeParse:
    """Returns queued ProposedClaim objects from messages.parse, recording each call."""

    def __init__(self, outputs):
        self.outputs = list(outputs)
        self.calls = []
        self.messages = SimpleNamespace(parse=self._parse)

    def _parse(self, **kw):
        self.calls.append(kw)
        return SimpleNamespace(parsed_output=self.outputs.pop(0), stop_reason="end_turn")


class FakeCreate:
    """Returns queued text responses from messages.create."""

    def __init__(self, texts):
        self.texts = list(texts)
        self.calls = []
        self.messages = SimpleNamespace(create=self._create)

    def _create(self, **kw):
        self.calls.append(kw)
        text = self.texts.pop(0)
        return SimpleNamespace(
            content=[SimpleNamespace(type="text", text=text)], stop_reason="end_turn"
        )


def _good_proposal(truth):
    return ProposedClaim(
        claim_type="differential_expression",
        name="DE_UP_0 is up in Mono under stim",
        condition_column="condition",
        sample_columns=["donor", "condition"],
        labels_column="cell_type",
        markers=mono_markers(truth)[:3],
        label_value="Mono",
        group_a="stim",
        group_b="ctrl",
        direction="up_in_a",
        genes=[truth["de_genes"][0]],
        rationale="donor x condition is the replicate; Mono has labels and markers",
    )


def test_describe_adata_shows_columns_values_and_matching_genes(sim, truth):
    text = describe_adata(
        sim, f"{truth['de_genes'][0]} is up in Mono cells after stim, not NOTAGENE"
    )
    assert "cells: 1,800" in text and "raw counts found in: .X" in text
    assert "condition: 2 values: ctrl" in text or "condition: 2 values: stim" in text
    assert "donor: 6 values" in text
    assert truth["de_genes"][0] in text.split("gene symbols from the claim")[1]
    assert "NOTAGENE" not in text.split("gene symbols from the claim")[1]


def test_parse_claim_returns_validated_spec(sim, truth):
    client = FakeParse([_good_proposal(truth)])
    res = parse_claim("DE_UP_0 up in Mono after stim", sim, client=client, model="fake-model")
    assert res.attempts == 1 and res.model == "fake-model"
    assert res.spec.claim_type == "differential_expression"
    assert res.spec.columns.sample_key == ["donor", "condition"]
    assert res.spec.genes == [truth["de_genes"][0]]
    sent = client.calls[0]
    assert sent["output_format"] is ProposedClaim
    assert "Dataset:" in sent["messages"][0]["content"]
    assert "cell_type" in sent["messages"][0]["content"]


def test_parse_claim_feeds_back_data_errors_once(sim, truth):
    bad = _good_proposal(truth).model_copy(
        update={"condition_column": "treatment", "genes": ["NOPE"]}
    )
    client = FakeParse([bad, _good_proposal(truth)])
    res = parse_claim("x", sim, client=client, model="m")
    assert res.attempts == 2 and len(client.calls) == 2
    retry_msgs = client.calls[1]["messages"]
    assert retry_msgs[1]["role"] == "assistant" and retry_msgs[2]["role"] == "user"
    assert "treatment" in retry_msgs[2]["content"] and "NOPE" in retry_msgs[2]["content"]


def test_parse_claim_gives_up_after_max_attempts(sim, truth):
    bad = _good_proposal(truth).model_copy(update={"group_a": "treated"})
    client = FakeParse([bad, bad])
    with pytest.raises(SpecError, match="after 2 attempts"):
        parse_claim("x", sim, client=client, model="m")


def test_proposed_claim_to_spec_drops_empty_fields():
    p = ProposedClaim(
        claim_type="population_exists", name="NK exist", markers=["GNLY"], rationale="r"
    )
    d = p.to_spec_dict()
    assert d["columns"] == {} and "comparison" not in d and "genes" not in d


# ---------------------------------------------------------------- summaries

PROV = {
    "claim": {"name": "ISG15 up in monocytes", "claim_type": "differential_expression"},
    "thresholds": {"robust_min_survival": 0.8, "fragile_min_survival": 0.2, "alpha": 0.05},
    "verdict": {
        "label": "FRAGILE",
        "n_runs": 21,
        "n_holds": 16,
        "survival": 16 / 21,
        "baseline_holds": True,
        "cell_level_survival": 1.0,
        "pseudoreplication_flag": False,
        "oat_flips": [{"params_changed": {"min_genes": 500}}],
        "breakers": [{"axis": "min_genes", "value": 500, "n_runs": 5, "n_holds": 0}],
    },
    "runs": [
        {
            "is_baseline": True,
            "match": {
                "n_cells": 5891,
                "marker_auroc": 0.953,
                "jaccard": 0.93,
                "n_samples_present": 16,
            },
            "details": {
                "per_gene": [
                    {
                        "gene": "ISG15",
                        "pb_log2fc": 7.26,
                        "pb_padj": 6.2e-61,
                        "cell_padj": 0.0,
                        "agreement": "both",
                    }
                ]
            },
        }
    ],
    "as_labelled": {"holds": True, "holds_cell_level": True},
}


def test_digest_contains_the_numbers_and_nothing_hidden():
    d = digest_from_provenance(PROV)
    assert "16 of 21" in d and "runs_failed: 5" in d and "min_genes=500" in d
    assert "FRAGILE" in d and "ISG15" in d and "7.26" in d


def test_number_guard_accepts_digest_numbers_and_rejects_others():
    d = digest_from_provenance(PROV)
    ok = (
        "The claim is FRAGILE: it held in 16 of 21 runs (76%) and min_genes=500 broke it "
        "in all 5 of those runs."
    )
    assert check_numbers(ok, d) == set()
    bad = "It failed in 6 runs and the log2 fold change was 7.3."
    assert check_numbers(bad, d) == {"6", "7.3"}


def test_summarize_retries_then_rejects():
    d_ok = "FRAGILE; held in 16 of 21 runs; min_genes=500 broke it."
    client = FakeCreate(["It failed 6 times.", d_ok])
    text, model_id = summarize_provenance(PROV, client=client, model="m")
    assert text == d_ok and model_id == "m" and len(client.calls) == 2
    assert "not in the digest" in client.calls[1]["messages"][-1]["content"]
    client2 = FakeCreate(["It failed 6 times.", "It failed 7 times."])
    with pytest.raises(SummaryRejected):
        summarize_provenance(PROV, client=client2, model="m")


def test_numbers_in_normalises():
    # numbers glued to letters (gene names such as ISG15, DE_UP_0) are not numbers
    got = numbers_in("16 of 21 (76%) p=0.050 and 1.0, but ISG15 and log2")
    assert got == {"16", "21", "76", "0.05", "1"}


# ---------------------------------------------------------------- no key


def test_parse_cli_without_key_points_to_hand_written_yaml(sim_path, tmp_path, monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("ANTHROPIC_AUTH_TOKEN", raising=False)
    monkeypatch.chdir(tmp_path)  # no .env here
    r = CliRunner().invoke(app, ["parse", "Mono exists", "--data", str(sim_path)])
    assert r.exit_code == 2
    assert "ANTHROPIC_API_KEY is not set" in r.output and "write the claim YAML by hand" in r.output


def test_run_with_summarize_but_no_key_still_completes(sim_path, tmp_path, monkeypatch, truth):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("ANTHROPIC_AUTH_TOKEN", raising=False)
    monkeypatch.chdir(tmp_path)
    claim = tmp_path / "c.yaml"
    claim.write_text(
        "claim_type: population_exists\nname: smoke\ncolumns: {sample: [donor, condition]}\n"
        f"target: {{markers: {json.dumps(mono_markers(truth)[:3])}}}\n"
    )
    r = CliRunner().invoke(
        app,
        [
            "run",
            "--data",
            str(sim_path),
            "--claim",
            str(claim),
            "--out",
            str(tmp_path),
            "--max-runs",
            "1",
            "-q",
            "--summarize",
        ],
    )
    assert r.exit_code == 0, r.output
    assert "no LLM summary" in r.output and (tmp_path / "smoke" / "report.html").exists()


def test_load_dotenv_only_sets_anthropic_keys(tmp_path, monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("OTHER", raising=False)
    env = tmp_path / ".env"
    env.write_text("# comment\nANTHROPIC_API_KEY='sk-test-not-real'\nOTHER=1\n")
    load_dotenv(env)
    import os

    assert os.environ["ANTHROPIC_API_KEY"] == "sk-test-not-real"
    assert "OTHER" not in os.environ
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
