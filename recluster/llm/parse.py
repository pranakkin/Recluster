"""Plain-English claim -> claim spec YAML, proposed by an LLM and verified against the data.

The model never sees the counts. It sees the claim text, the ``.obs`` columns with their
values, and which gene symbols from the claim exist in ``var_names``. Its proposal is
converted to a :class:`ClaimSpec`, validated, and checked against the AnnData; validation
errors are fed back once. The user confirms the YAML before anything runs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Literal

from anndata import AnnData
from pydantic import BaseModel, Field, ValidationError

from recluster.data import extract_counts
from recluster.llm.common import LLMUnavailable, get_client, model_name
from recluster.spec import ClaimSpec, SpecError

MAX_VALUES = 20


class ProposedClaim(BaseModel):
    """Flat, schema-friendly shape the model fills in. Empty string means "not applicable"."""

    claim_type: Literal["population_exists", "composition_shift", "differential_expression"]
    name: str = Field(description="one short sentence stating the claim")
    condition_column: str = Field(default="", description=".obs column holding the condition")
    sample_columns: list[str] = Field(
        default_factory=list,
        description="obs column(s) identifying the biological replicate; several are combined",
    )
    batch_column: str = Field(default="", description=".obs column for technical batch, if any")
    labels_column: str = Field(default="", description=".obs column with existing cell-type labels")
    markers: list[str] = Field(
        default_factory=list, description="marker genes defining the population"
    )
    negative_markers: list[str] = Field(default_factory=list)
    label_value: str = Field(default="", description="value in labels_column naming the population")
    group_a: str = Field(default="", description="condition value for group A")
    group_b: str = Field(default="", description="condition value for group B")
    direction: str = Field(
        default="",
        description=(
            "composition: higher_in_a|lower_in_a|different; expression: up_in_a|down_in_a|either"
        ),
    )
    genes: list[str] = Field(
        default_factory=list, description="genes the claim is about (expression only)"
    )
    require: Literal["all", "any"] = "all"
    rationale: str = Field(description="one sentence on how the columns and values were chosen")

    def to_spec_dict(self) -> dict[str, Any]:
        cols: dict[str, Any] = {}
        if self.condition_column:
            cols["condition"] = self.condition_column
        if self.sample_columns:
            cols["sample"] = (
                self.sample_columns if len(self.sample_columns) > 1 else self.sample_columns[0]
            )
        if self.batch_column:
            cols["batch"] = self.batch_column
        if self.labels_column:
            cols["labels"] = self.labels_column
        target: dict[str, Any] = {}
        if self.markers:
            target["markers"] = self.markers
        if self.negative_markers:
            target["negative_markers"] = self.negative_markers
        if self.label_value:
            target["label_value"] = self.label_value
        d: dict[str, Any] = {
            "claim_type": self.claim_type,
            "name": self.name,
            "columns": cols,
            "target": target,
        }
        if self.claim_type != "population_exists":
            d["comparison"] = {
                "group_a": self.group_a,
                "group_b": self.group_b,
                "direction": self.direction,
            }
        if self.claim_type == "differential_expression":
            d["genes"] = self.genes
            d["require"] = self.require
        return d


@dataclass
class ParseResult:
    spec: ClaimSpec
    rationale: str
    model: str
    attempts: int
    data_summary: str


# ----------------------------------------------------------------- data view


def describe_adata(adata: AnnData, claim_text: str = "") -> str:
    """What the model is allowed to know about the data: columns, values, matching genes."""
    lines = [f"cells: {adata.n_obs:,}; genes: {adata.n_vars:,}"]
    try:
        _, source = extract_counts(adata)
        lines.append(f"raw counts found in: {source}")
    except Exception as exc:  # noqa: BLE001 - informational only
        lines.append(f"raw counts: not found ({exc.__class__.__name__})")
    lines.append("obs columns:")
    for col in adata.obs.columns:
        s = adata.obs[col]
        n_unique = int(s.nunique())
        if s.dtype.name in ("category", "object", "bool", "string") or n_unique <= 30:
            counts = s.astype(str).value_counts()
            shown = ", ".join(f"{v} ({n})" for v, n in counts.head(MAX_VALUES).items())
            more = f", ... {n_unique - MAX_VALUES} more" if n_unique > MAX_VALUES else ""
            lines.append(f"  - {col}: {n_unique} values: {shown}{more}")
        else:
            lines.append(f"  - {col}: numeric, {n_unique} distinct values")
    if claim_text:
        tokens = sorted({t for t in re.findall(r"[A-Za-z][A-Za-z0-9._\-]{1,}", claim_text)})
        names = set(adata.var_names.astype(str))
        if adata.raw is not None:
            names |= set(adata.raw.var_names.astype(str))
        found = [t for t in tokens if t in names]
        lines.append(
            "gene symbols from the claim that exist in var_names: "
            + (", ".join(found) if found else "none")
        )
        lines.append("any other gene you name must be spelled exactly as in var_names")
    return "\n".join(lines)


SYSTEM = """You turn a plain-English claim about a single-cell RNA-seq dataset into a structured \
claim specification for Recluster, a tool that stress-tests such claims. You will be given the \
claim and a description of the dataset's .obs columns and values. Rules:
- Use column names and values exactly as listed. Never invent a column or a value.
- claim_type: population_exists (a population exists), composition_shift (a population's \
share differs between two conditions), differential_expression (named genes change in a \
population between conditions).
- sample_columns must identify the biological replicate (donor, patient, mouse, sample). If \
donors appear in both conditions, list both the donor column and the condition column so each \
replicate is nested in one condition. population_exists claims may also set sample_columns \
when a replicate column exists.
- Define the population by marker genes (markers, optional negative_markers) and/or by an \
existing label (labels_column + label_value) when the dataset has cell-type labels that match \
the claim.
- group_a is the condition the claim says is higher/up; group_b the comparison. Directions: \
composition higher_in_a / lower_in_a / different; expression up_in_a / down_in_a / either.
- genes only for differential_expression; require is "all" unless the claim says "any".
- name: one short sentence restating the claim. rationale: one sentence on your choices.
- If the claim cannot be expressed (no replicate column for a comparison claim, no condition \
column), still fill the closest spec and explain the gap in rationale."""


def parse_claim(
    claim_text: str,
    adata: AnnData,
    *,
    client: Any = None,
    model: str | None = None,
    max_attempts: int = 2,
    audit: bool = False,
) -> ParseResult:
    """Propose a spec, validate it against the data, feed errors back once."""
    client = client or get_client()
    model_id = model_name(model)
    summary = describe_adata(adata, claim_text)
    user = f"Claim: {claim_text}\n\nDataset:\n{summary}"
    messages: list[dict[str, Any]] = [{"role": "user", "content": user}]
    last_error = ""
    for attempt in range(1, max_attempts + 1):
        response = client.messages.parse(
            model=model_id,
            max_tokens=4000,
            system=SYSTEM,
            messages=messages,
            output_format=ProposedClaim,
        )
        if getattr(response, "stop_reason", "") == "refusal":
            raise LLMUnavailable("the model declined to parse this claim; write the YAML by hand")
        proposed: ProposedClaim = response.parsed_output
        try:
            spec = ClaimSpec.model_validate(proposed.to_spec_dict())
            spec.check_against(adata, audit=audit)
            return ParseResult(spec, proposed.rationale, model_id, attempt, summary)
        except (ValidationError, SpecError) as exc:
            last_error = str(exc)
            messages += [
                {"role": "assistant", "content": proposed.model_dump_json()},
                {
                    "role": "user",
                    "content": (
                        "That specification was rejected when checked against the data:\n"
                        f"{last_error}\nFix it using only the listed columns, values and genes."
                    ),
                },
            ]
    raise SpecError(
        f"the proposed spec still does not match the data after {max_attempts} attempts:\n"
        f"{last_error}"
    )


def yaml_with_header(result: ParseResult, claim_text: str) -> str:
    """The YAML to show and save, labelled as LLM output."""
    header = (
        f"# Proposed by an LLM ({result.model}) from the plain-English claim:\n"
        f"#   {claim_text}\n"
        f"# Rationale (LLM): {result.rationale}\n"
        "# Checked against the data: all columns, values and genes exist. Review before running.\n"
    )
    return header + result.spec.to_yaml()
