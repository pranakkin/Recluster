"""Claim spec validation: internal consistency and agreement with the data."""

from __future__ import annotations

import pytest
import yaml
from pydantic import ValidationError

from recluster.spec import ClaimSpec, SpecError

BASE_COLS = {"condition": "condition", "sample": ["donor", "condition"], "labels": "cell_type"}


def _de(**over) -> dict:
    d = {
        "claim_type": "differential_expression",
        "name": "x",
        "columns": BASE_COLS,
        "target": {"label_value": "Mono"},
        "comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "up_in_a"},
        "genes": ["DE_UP_0"],
    }
    d.update(over)
    return d


def test_valid_specs_round_trip(tmp_path):
    spec = ClaimSpec.model_validate(_de())
    p = tmp_path / "c.yaml"
    p.write_text(spec.to_yaml())
    again = ClaimSpec.from_yaml(p)
    assert again == spec
    assert yaml.safe_load(spec.to_yaml())["claim_type"] == "differential_expression"


def test_population_spec_minimal():
    spec = ClaimSpec.model_validate(
        {"claim_type": "population_exists", "name": "x", "target": {"markers": ["CD14"]}}
    )
    assert spec.columns.sample_key is None
    assert spec.comparison is None


@pytest.mark.parametrize(
    "bad, msg",
    [
        ({"target": {"min_cells": 5}}, "markers"),
        ({"comparison": None}, "requires `comparison`"),
        ({"genes": []}, "at least one gene"),
        (
            {"comparison": {"group_a": "stim", "group_b": "ctrl", "direction": "higher_in_a"}},
            "direction",
        ),
        ({"comparison": {"group_a": "stim", "group_b": "stim", "direction": "up_in_a"}}, "differ"),
        ({"columns": {"condition": "condition", "labels": "cell_type"}}, "columns.sample"),
        ({"columns": {"condition": "condition", "sample": "s"}}, "requires `columns.labels`"),
        ({"require": "most"}, "require"),
        ({"unknown_field": 1}, "unknown_field"),
    ],
)
def test_inconsistent_specs_rejected(bad, msg):
    with pytest.raises(ValidationError) as exc:
        ClaimSpec.model_validate(_de(**bad))
    assert msg in str(exc.value)


def test_population_rejects_comparison_and_genes():
    with pytest.raises(ValidationError):
        ClaimSpec.model_validate(
            {
                "claim_type": "population_exists",
                "name": "x",
                "target": {"markers": ["A"]},
                "comparison": {"group_a": "a", "group_b": "b", "direction": "different"},
            }
        )
    with pytest.raises(ValidationError):
        ClaimSpec.model_validate(
            {
                "claim_type": "population_exists",
                "name": "x",
                "target": {"markers": ["A"]},
                "genes": ["G"],
            }
        )


def test_check_against_passes_on_matching_data(sim):
    ClaimSpec.model_validate(_de()).check_against(sim)


def test_check_against_reports_every_problem(sim):
    spec = ClaimSpec.model_validate(
        _de(
            columns={"condition": "nope", "sample": "donor", "labels": "cell_type", "batch": "zzz"},
            target={"label_value": "Granulocyte"},
            genes=["NOT_A_GENE"],
        )
    )
    with pytest.raises(SpecError) as exc:
        spec.check_against(sim)
    text = str(exc.value)
    assert "columns.condition='nope'" in text
    assert "columns.batch='zzz'" in text
    assert "Granulocyte" in text
    assert "NOT_A_GENE" in text


def test_check_against_requires_samples_nested_in_condition(sim):
    spec = ClaimSpec.model_validate(_de(columns={**BASE_COLS, "sample": "donor"}))
    with pytest.raises(SpecError, match="span more than one condition"):
        spec.check_against(sim)


def test_check_against_unknown_group(sim):
    spec = ClaimSpec.model_validate(
        _de(comparison={"group_a": "treated", "group_b": "ctrl", "direction": "up_in_a"})
    )
    with pytest.raises(SpecError, match="'treated' not in"):
        spec.check_against(sim)


def test_audit_requires_labels(sim):
    spec = ClaimSpec.model_validate(
        {"claim_type": "population_exists", "name": "x", "target": {"markers": ["Mono_marker_0"]}}
    )
    with pytest.raises(SpecError, match="audit mode requires"):
        spec.check_against(sim, audit=True)
