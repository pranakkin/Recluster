"""Claim specification: pydantic models, YAML loading, and validation against an AnnData.

A claim spec says *what* is being asserted and *where* in the data to look. It never
contains results. The three claim types are:

- ``population_exists``: a population defined by markers and/or an existing label exists
  as a distinct cluster.
- ``composition_shift``: the proportion of a population differs between two conditions.
- ``differential_expression``: named genes change in a population between two conditions.
"""

from __future__ import annotations

import difflib
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

import pandas as pd
import yaml
from anndata import AnnData
from pydantic import BaseModel, ConfigDict, Field, model_validator

ClaimType = Literal["population_exists", "composition_shift", "differential_expression"]

CompositionDirection = Literal["higher_in_a", "lower_in_a", "different"]
DEDirection = Literal["up_in_a", "down_in_a", "either"]


class SpecError(ValueError):
    """Raised when a claim spec is internally inconsistent or does not match the data."""


class Columns(BaseModel):
    """Names of the relevant ``.obs`` columns."""

    model_config = ConfigDict(extra="forbid")

    condition: str | None = None
    sample: str | list[str] | None = None
    batch: str | None = None
    labels: str | None = None
    collection_time: str | None = Field(
        default=None,
        description="column with the collection time of each sample (hours 0-24 or 'HH:MM')",
    )

    @property
    def sample_key(self) -> list[str] | None:
        """Sample column(s) as a list, or None."""
        if self.sample is None:
            return None
        return [self.sample] if isinstance(self.sample, str) else list(self.sample)


class Target(BaseModel):
    """How to find the population a claim is about."""

    model_config = ConfigDict(extra="forbid")

    markers: list[str] = Field(default_factory=list)
    negative_markers: list[str] = Field(default_factory=list)
    marker_logic: Literal["any", "all"] = Field(
        default="any",
        description=(
            "'any': markers are alternative evidence for one type (mean score); "
            "'all': the population co-expresses every marker (minimum score)"
        ),
    )
    label_value: str | None = None
    min_cells: int = Field(default=20, ge=1)

    @model_validator(mode="after")
    def _at_least_one_definition(self) -> Target:
        if not self.markers and self.label_value is None:
            raise ValueError("target needs `markers` and/or `label_value`")
        overlap = set(self.markers) & set(self.negative_markers)
        if overlap:
            raise ValueError(f"genes listed as both marker and negative marker: {sorted(overlap)}")
        return self


class Comparison(BaseModel):
    """Two condition values and the asserted direction of the effect."""

    model_config = ConfigDict(extra="forbid")

    group_a: str
    group_b: str
    direction: str

    @model_validator(mode="after")
    def _distinct_groups(self) -> Comparison:
        if self.group_a == self.group_b:
            raise ValueError("group_a and group_b must differ")
        return self


class ClaimSpec(BaseModel):
    """A structured, machine-checkable claim about a single-cell dataset."""

    model_config = ConfigDict(extra="forbid")

    claim_type: ClaimType
    name: str = Field(min_length=1)
    columns: Columns = Field(default_factory=Columns)
    target: Target
    comparison: Comparison | None = None
    genes: list[str] = Field(default_factory=list)
    require: Literal["all", "any"] = "all"
    tissue: str | None = Field(
        default=None,
        description="tissue name for the public day-night gene list (e.g. 'Whole Blood')",
    )

    @model_validator(mode="after")
    def _consistent(self) -> ClaimSpec:
        needs_comparison = self.claim_type in ("composition_shift", "differential_expression")
        if needs_comparison:
            if self.comparison is None:
                raise ValueError(f"{self.claim_type} requires `comparison`")
            if self.columns.condition is None:
                raise ValueError(f"{self.claim_type} requires `columns.condition`")
            if self.columns.sample is None:
                raise ValueError(
                    f"{self.claim_type} requires `columns.sample` (biological replicates)"
                )
        if self.claim_type == "population_exists" and self.comparison is not None:
            raise ValueError("population_exists does not take `comparison`")
        if self.claim_type == "differential_expression":
            if not self.genes:
                raise ValueError("differential_expression requires at least one gene")
            allowed = DEDirection.__args__  # type: ignore[attr-defined]
        elif self.claim_type == "composition_shift":
            allowed = CompositionDirection.__args__  # type: ignore[attr-defined]
        else:
            allowed = ()
        if self.comparison is not None and self.comparison.direction not in allowed:
            raise ValueError(
                f"direction for {self.claim_type} must be one of {list(allowed)}, "
                f"got {self.comparison.direction!r}"
            )
        if self.claim_type != "differential_expression" and self.genes:
            raise ValueError("`genes` is only used by differential_expression")
        if self.target.label_value is not None and self.columns.labels is None:
            raise ValueError("target.label_value requires `columns.labels`")
        return self

    # ---------------------------------------------------------------- IO

    @classmethod
    def from_yaml(cls, path: str | Path) -> ClaimSpec:
        """Load a spec from a YAML file."""
        with open(path) as fh:
            data = yaml.safe_load(fh)
        if data is None:
            raise SpecError(f"{path} is empty")
        if not isinstance(data, dict):
            raise SpecError(
                f"{path}: expected `key: value` lines at the top level "
                "(start from examples/claims/)"
            )
        return cls.model_validate(data)

    def to_yaml(self) -> str:
        """Serialise to YAML (used by the LLM parser and for provenance)."""
        return yaml.safe_dump(self.model_dump(exclude_none=True), sort_keys=False)

    # ---------------------------------------------------------- data check

    def check_against(self, adata: AnnData, *, audit: bool = False) -> None:
        """Verify every named column, value and gene exists in ``adata``.

        Raises :class:`SpecError` listing all problems found.
        """
        problems: list[str] = []
        obs = adata.obs
        cols = self.columns
        if audit and cols.labels is None:
            problems.append(
                "audit mode requires `columns.labels` "
                "(the .obs column holding the other tool's cell labels)"
            )

        for role, col in (
            ("condition", cols.condition),
            ("batch", cols.batch),
            ("labels", cols.labels),
        ):
            if col is not None and col not in obs.columns:
                problems.append(_missing_column(f"columns.{role}", col, obs.columns))
        for col in cols.sample_key or []:
            if col not in obs.columns:
                problems.append(_missing_column("columns.sample", col, obs.columns))
        if cols.collection_time is not None:
            if cols.collection_time not in obs.columns:
                problems.append(
                    _missing_column("columns.collection_time", cols.collection_time, obs.columns)
                )
            else:
                from recluster.hours import parse_hours

                bad = parse_hours(obs[cols.collection_time]).isna().mean()
                if bad > 0:
                    problems.append(
                        f"{bad:.0%} of collection_time values could not be read as hours or HH:MM"
                    )

        if (
            cols.labels is not None
            and cols.labels in obs.columns
            and self.target.label_value is not None
        ):
            values = obs[cols.labels].astype(str)
            if self.target.label_value not in set(values):
                problems.append(
                    f"target.label_value={self.target.label_value!r} not found in "
                    f".obs[{cols.labels!r}]{_near(self.target.label_value, values.unique())} "
                    f"(have: {_preview(values.unique())})"
                )

        if self.comparison is not None and cols.condition in obs.columns:
            present = set(obs[cols.condition].astype(str))
            for grp in (self.comparison.group_a, self.comparison.group_b):
                if grp not in present:
                    problems.append(
                        f"comparison group {grp!r} not in .obs[{cols.condition!r}]"
                        f"{_near(grp, sorted(present))} (have: {_preview(sorted(present))})"
                    )
            if (
                cols.condition is not None
                and cols.sample_key
                and all(c in obs.columns for c in cols.sample_key)
            ):
                problems.extend(_check_nesting(obs, cols.sample_key, cols.condition))

        gene_names = set(adata.var_names)
        if adata.raw is not None:
            gene_names |= set(adata.raw.var_names)
        any_missing = False
        for role, genes in (
            ("target.markers", self.target.markers),
            ("target.negative_markers", self.target.negative_markers),
            ("genes", self.genes),
        ):
            missing = [g for g in genes if g not in gene_names]
            if missing:
                any_missing = True
                problems.append(
                    f"{role}: {len(missing)} gene(s) not found in the data's gene names "
                    f"(var_names): {missing}{_gene_hints(missing, gene_names)}"
                )
        if any_missing and _mostly_ensembl(gene_names):
            problems.append(
                "most gene names in this file are Ensembl IDs (ENS...): name genes by their "
                "Ensembl ID in the claim, or set gene symbols as var_names"
            )

        if problems:
            raise SpecError(
                "claim spec does not match the data (fix these in the YAML and check again):"
                "\n  - " + "\n  - ".join(problems)
            )


def _preview(values: object, n: int = 8) -> str:
    vals = list(values)  # type: ignore[call-overload]
    shown = ", ".join(repr(v) for v in vals[:n])
    return shown + (", ..." if len(vals) > n else "")


def _closest(value: str, options: Iterable[object], cutoff: float = 0.6) -> str | None:
    """The option closest to ``value``: a case-insensitive exact match first, then difflib."""
    opts = [str(o) for o in options]
    folded = [o for o in opts if o.lower() == value.lower()]
    if folded:
        return folded[0]
    close = difflib.get_close_matches(value, opts, n=1, cutoff=cutoff)
    if close:
        return close[0]
    # an abbreviation or a longer spelling, e.g. 'stimulated' for 'stim'
    low = value.lower()
    prefixed = [
        o
        for o in opts
        if min(len(o), len(value)) >= 3 and (low.startswith(o.lower()) or o.lower().startswith(low))
    ]
    return prefixed[0] if len(prefixed) == 1 else None


def _near(value: str, options: Iterable[object]) -> str:
    """'. Did you mean ...?' for the closest option, or '' when nothing is close."""
    best = _closest(value, options)
    if best is None or best == value:
        return ""
    case = " (values are case-sensitive)" if best.lower() == value.lower() else ""
    return f". Did you mean {best!r}?{case}"


def _missing_column(field: str, col: str, columns: Iterable[object]) -> str:
    cols = list(columns)
    return f"{field}={col!r} not in .obs{_near(col, cols)} (have: {_preview(cols)})"


def _gene_hints(missing: list[str], gene_names: set[str], max_genes: int = 10) -> str:
    """Suggestions for misspelt gene names: same name in another case, or a close spelling."""
    by_lower: dict[str, str] = {}
    for g in gene_names:
        by_lower.setdefault(str(g).lower(), str(g))
    pool = sorted(str(g) for g in gene_names)
    hints = []
    for g in missing[:max_genes]:
        best = by_lower.get(g.lower())
        if best is None:
            close = difflib.get_close_matches(g, pool, n=1, cutoff=0.85)
            best = close[0] if close else None
        if best is not None:
            hints.append(f"{g!r} -> {best!r}")
    if not hints:
        return ""
    return ". Did you mean " + ", ".join(hints) + "? (gene names are case-sensitive)"


def _mostly_ensembl(gene_names: set[str]) -> bool:
    if not gene_names:
        return False
    n_ens = sum(1 for g in gene_names if str(g).startswith("ENS"))
    return n_ens > len(gene_names) / 2


def _check_nesting(obs: pd.DataFrame, sample_key: list[str], condition: str) -> list[str]:
    """Each sample must belong to exactly one condition."""
    key = make_sample_key(obs, sample_key)
    n_conditions = obs.groupby(key, observed=True)[condition].nunique()
    bad = n_conditions[n_conditions > 1]
    if len(bad) == 0:
        return []
    return [
        f"{len(bad)} sample(s) span more than one condition (e.g. {_preview(bad.index, 3)}). "
        f"Replicates must be nested in condition; use a combined key such as "
        f"columns.sample: [{', '.join(sample_key)}, {condition}]"
    ]


def make_sample_key(obs: pd.DataFrame, sample_key: list[str]) -> pd.Series:
    """Combine one or more obs columns into a single string replicate id."""
    parts = [obs[c].astype(str) for c in sample_key]
    out = parts[0]
    for p in parts[1:]:
        out = out + "|" + p
    return out.rename("sample")
