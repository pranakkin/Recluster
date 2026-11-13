"""Command-line interface."""

from __future__ import annotations

import difflib
import re
import sys
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Any

import typer

from recluster import __version__

if TYPE_CHECKING:
    from pydantic import BaseModel, ValidationError

    from recluster.config import GridConfig, Thresholds
    from recluster.spec import ClaimSpec

EPILOG = (
    "**Start here:**\n\n"
    "`recluster demo` runs two example claims on a public dataset\n\n"
    "copy examples/claims/TEMPLATE.yaml and edit it for your data\n\n"
    "`recluster check --data X.h5ad --claim claim.yaml` checks the claim against the data\n\n"
    "`recluster run --data X.h5ad --claim claim.yaml` runs the stress test"
)

app = typer.Typer(
    name="recluster",
    help=(
        "Stress-test a single-cell RNA-seq claim across a grid of defensible analysis choices "
        "and write an HTML report with a verdict: ROBUST, FRAGILE or NOT SUPPORTED."
    ),
    epilog=EPILOG,
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode="markdown",
)

DataOpt = Annotated[
    Path,
    typer.Option(
        "--data",
        help="AnnData .h5ad file with raw integer counts (in layers['counts'], .raw.X or .X)",
        exists=True,
        dir_okay=False,
    ),
]
ClaimOpt = Annotated[
    Path,
    typer.Option(
        "--claim",
        help="claim YAML (copy one from examples/claims/)",
        exists=True,
        dir_okay=False,
    ),
]
OutOpt = Annotated[
    Path,
    typer.Option(
        "--out",
        help="results folder; each claim gets a subfolder named after it with report.html, "
        "provenance.json, runs.csv and claim.yaml",
    ),
]
GridOpt = Annotated[
    Path | None,
    typer.Option(
        "--grid",
        help="YAML overriding grid axes and settings; keys left out keep their defaults; "
        "the report prints the values used",
        exists=True,
        dir_okay=False,
    ),
]
ThrOpt = Annotated[
    Path | None,
    typer.Option(
        "--thresholds",
        help="YAML overriding verdict thresholds; keys left out keep their defaults; "
        "the report prints the values used",
        exists=True,
        dir_okay=False,
    ),
]
MaxRunsOpt = Annotated[
    int | None,
    typer.Option(
        "--max-runs",
        help="run at most this many grid points, for a smoke test only: the verdict from a "
        "capped grid is not meaningful [default: the whole grid]",
    ),
]
RandomOpt = Annotated[
    int | None,
    typer.Option(
        "--n-random",
        help="random combinations added to the one-at-a-time runs "
        "[default: 6, or the --grid file's value]",
    ),
]
QuietOpt = Annotated[bool, typer.Option("--quiet", "-q", help="hide per-run progress")]
JobsOpt = Annotated[
    int,
    typer.Option(
        "--jobs",
        "-j",
        min=1,
        help="worker processes for the grid runs; results are identical to --jobs 1, "
        "only faster (each worker holds a copy of the counts)",
    ),
]
SummarizeOpt = Annotated[
    bool,
    typer.Option(
        "--summarize",
        help="add an LLM-written summary of the computed results (optional; needs "
        "`uv sync --extra llm` and ANTHROPIC_API_KEY)",
    ),
]

EXAMPLES_HINT = "Examples to copy: examples/claims/"


# ------------------------------------------------------------------ readable file errors


def _loc(parts: tuple[Any, ...]) -> str:
    out = ""
    for p in parts:
        out += f"[{p}]" if isinstance(p, int) else (f".{p}" if out else str(p))
    return out


def _describe_error(err: dict[str, Any], models: dict[str, type[BaseModel]]) -> str:
    """One plain sentence for one pydantic error, with a nearest-match hint where possible."""
    loc = tuple(err.get("loc", ()))
    where = _loc(loc)
    kind = err.get("type", "")
    if kind == "extra_forbidden" and loc:
        key = str(loc[-1])
        parent = ".".join(str(p) for p in loc[:-1] if not isinstance(p, int))
        model = models.get(parent)
        allowed = list(model.model_fields) if model is not None else []
        close = difflib.get_close_matches(key, allowed, n=1, cutoff=0.6)
        hint = f" Did you mean {close[0]!r}?" if close else ""
        if not close:
            # the key may belong one level down, e.g. `resolution` under `axes`
            elsewhere = [
                f"{prefix}.{m[0]}"
                for prefix, other in models.items()
                if prefix and prefix != parent
                for m in [difflib.get_close_matches(key, list(other.model_fields), 1, 0.8)]
                if m
            ]
            if elsewhere:
                hint = " Did you mean " + " or ".join(repr(e) for e in elsewhere[:2]) + "?"
        listed = f" (allowed here: {', '.join(allowed)})" if allowed else ""
        return f"unknown key `{where}`.{hint}{listed}"
    if kind == "literal_error":
        allowed = re.findall(r"'([^']*)'", str(err.get("ctx", {}).get("expected", "")))
        value = err.get("input")
        close = difflib.get_close_matches(str(value), allowed, n=1, cutoff=0.6)
        hint = f" Did you mean {close[0]!r}?" if close else ""
        return f"`{where}`: {value!r} is not allowed.{hint} (allowed: {', '.join(allowed)})"
    if kind == "missing":
        return f"`{where}` is required"
    msg = str(err.get("msg", "")).removeprefix("Value error, ")
    value = err.get("input")
    got = f" (got {value!r})" if isinstance(value, str | int | float | bool) else ""
    return f"`{where}`: {msg}{got}" if where else msg


def _format_validation(exc: ValidationError, models: dict[str, type[BaseModel]]) -> list[str]:
    return [_describe_error(dict(e), models) for e in exc.errors()]


def _read_yaml_model(path: Path, header: str, load: Any, models: dict[str, Any]) -> Any:
    """Call ``load(path)``; on failure print readable bullets under ``header`` and exit 2."""
    import yaml
    from pydantic import ValidationError

    from recluster.spec import SpecError

    problems: list[str]
    if path.suffix.lower() in (".h5ad", ".h5"):
        problems = [
            f"{path.name} looks like a data file: --claim takes the claim YAML and --data the "
            ".h5ad (are the two swapped?)"
        ]
    else:
        try:
            return load(path)
        except yaml.YAMLError as exc:
            mark = getattr(exc, "problem_mark", None)
            problem = getattr(exc, "problem", None) or str(exc)
            at = f" at line {mark.line + 1}, column {mark.column + 1}" if mark is not None else ""
            problems = [
                f"not valid YAML{at}: {problem}. Check for an unclosed quote or bracket, "
                "or tabs used for indentation."
            ]
        except UnicodeDecodeError:
            problems = ["not a text file; expected a YAML file"]
        except SpecError as exc:
            problems = [str(exc)]
        except ValidationError as exc:
            problems = _format_validation(exc, models)
    lines = [f"{header} {path}:", *(f"  - {p}" for p in problems)]
    if header.startswith("invalid claim spec"):
        lines.append(EXAMPLES_HINT)
    typer.secho("\n".join(lines), err=True, fg="red")
    raise typer.Exit(2)


def _claim_models() -> dict[str, Any]:
    from recluster.spec import ClaimSpec, Columns, Comparison, Target

    return {"": ClaimSpec, "columns": Columns, "target": Target, "comparison": Comparison}


def _load_claim(path: Path) -> ClaimSpec:
    """Load a claim YAML, or print what is wrong with it (with suggestions) and exit 2."""
    from recluster.spec import ClaimSpec

    spec: ClaimSpec = _read_yaml_model(
        path, "invalid claim spec", ClaimSpec.from_yaml, _claim_models()
    )
    return spec


def _load_grid(path: Path | None) -> GridConfig:
    from recluster.config import GridAxes, GridConfig, PipelineParams

    if path is None:
        return GridConfig()
    models = {"": GridConfig, "axes": GridAxes, "baseline": PipelineParams}
    cfg: GridConfig = _read_yaml_model(path, "invalid grid file", GridConfig.from_yaml, models)
    return cfg


def _load_thresholds(path: Path | None) -> Thresholds:
    from recluster.config import Thresholds

    if path is None:
        return Thresholds()
    thr: Thresholds = _read_yaml_model(
        path, "invalid thresholds file", Thresholds.from_yaml, {"": Thresholds}
    )
    return thr


# ------------------------------------------------------------------ commands


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            help="print the version and exit",
            callback=_version_callback,
            is_eager=True,
        ),
    ] = False,
) -> None:
    pass


def _execute(
    data: Path,
    claim: Path,
    out: Path,
    grid: Path | None,
    thresholds: Path | None,
    max_runs: int | None,
    n_random: int | None,
    quiet: bool,
    audit: bool,
    summarize: bool = False,
    jobs: int = 1,
) -> None:
    from recluster.data import DataError
    from recluster.engine import run_claim, slugify
    from recluster.spec import SpecError

    spec = _load_claim(claim)
    cfg = _load_grid(grid)
    if max_runs is not None:
        cfg = cfg.model_copy(update={"max_runs": max_runs})
    if n_random is not None:
        cfg = cfg.model_copy(update={"n_random": n_random})
    thr = _load_thresholds(thresholds)
    try:
        result = run_claim(
            data,
            spec,
            out_dir=out,
            grid=cfg,
            thresholds=thr,
            audit=audit,
            jobs=jobs,
            progress=None if quiet else (lambda m: print(m, file=sys.stderr, flush=True)),
        )
    except (SpecError, DataError) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    v = result.verdict
    typer.echo(f"{v.label}\t{v.n_holds}/{v.n_runs} runs\t{spec.name}")
    typer.echo(f"report: {out / slugify(spec.name) / 'report.html'}")
    if summarize:
        _attach_summary(result, out, spec)


def _attach_summary(result: object, out: Path, spec: object) -> None:
    """Optional LLM summary after a run; failures never affect the computed outputs."""
    from recluster.engine import slugify
    from recluster.llm.common import LLMUnavailable
    from recluster.llm.summarize import SummaryRejected, summarize_provenance, summary_markdown
    from recluster.report.html import write_report

    folder = out / slugify(spec.name)  # type: ignore[attr-defined]
    try:
        text, model_id = summarize_provenance(result.provenance)  # type: ignore[attr-defined]
    except LLMUnavailable as exc:
        typer.secho(f"no LLM summary: {exc}", err=True, fg="yellow")
        return
    except SummaryRejected as exc:
        typer.secho(f"LLM summary discarded: {exc}", err=True, fg="yellow")
        return
    (folder / "llm_summary.md").write_text(summary_markdown(text, model_id))
    write_report(result, folder / "report.html", llm_summary=text)  # type: ignore[arg-type]
    typer.echo(f"\nLLM summary of computed results ({model_id}):\n{text}")


@app.command()
def check(
    data: DataOpt,
    claim: ClaimOpt,
    audit: Annotated[
        bool,
        typer.Option(
            "--audit",
            help="check for `recluster audit`: the data is a processed object with its own "
            "cell labels",
        ),
    ] = False,
) -> None:
    """Check a claim YAML against the data without running anything: columns, values, genes,
    raw counts.
    """
    from recluster.data import AUDIT_COUNTS_IN_X, DataError, extract_counts, load_adata
    from recluster.spec import SpecError

    spec = _load_claim(claim)
    try:
        adata = load_adata(data)
        spec.check_against(adata, audit=audit)
        counts, source = extract_counts(adata)
        if audit and source == ".X":
            raise DataError(AUDIT_COUNTS_IN_X)
    except (SpecError, DataError) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    command = "audit" if audit else "run"
    typer.echo(
        f"OK: {spec.name!r} matches {data} ({counts.n_obs} cells, {counts.n_vars} genes; "
        f"raw counts in {source})."
    )
    typer.echo(f"Next: recluster {command} --data {data} --claim {claim}")


@app.command()
def run(
    data: DataOpt,
    claim: ClaimOpt,
    out: OutOpt = Path("results"),
    grid: GridOpt = None,
    thresholds: ThrOpt = None,
    max_runs: MaxRunsOpt = None,
    n_random: RandomOpt = None,
    quiet: QuietOpt = False,
    summarize: SummarizeOpt = False,
    jobs: JobsOpt = 1,
) -> None:
    """Stress-test a claim: re-run the analysis from raw counts under each grid setting and
    give a verdict. Writes `OUT/<claim-name>/report.html`.
    """
    _execute(
        data,
        claim,
        out,
        grid,
        thresholds,
        max_runs,
        n_random,
        quiet,
        audit=False,
        summarize=summarize,
        jobs=jobs,
    )


@app.command()
def audit(
    data: DataOpt,
    claim: ClaimOpt,
    out: OutOpt = Path("results"),
    grid: GridOpt = None,
    thresholds: ThrOpt = None,
    max_runs: MaxRunsOpt = None,
    n_random: RandomOpt = None,
    quiet: QuietOpt = False,
    summarize: SummarizeOpt = False,
    jobs: JobsOpt = 1,
) -> None:
    """Check another tool's cell labels: re-derive from the processed object's stored raw
    counts and compare. Needs raw counts in layers['counts'] or .raw and `columns.labels` in
    the claim.
    """
    _execute(
        data,
        claim,
        out,
        grid,
        thresholds,
        max_runs,
        n_random,
        quiet,
        audit=True,
        summarize=summarize,
        jobs=jobs,
    )


class BenchPart(StrEnum):
    """Which part of the benchmark to run."""

    kang = "kang"
    planted = "planted"
    all = "all"


@app.command()
def bench(
    data_dir: Annotated[
        Path, typer.Option("--data-dir", help="where datasets are downloaded")
    ] = Path("data"),
    out: Annotated[Path, typer.Option("--out", help="where the benchmark reports go")] = Path(
        "bench/out"
    ),
    part: Annotated[BenchPart, typer.Option("--part", help="which part to run")] = BenchPart.all,
    jobs: Annotated[int, typer.Option("--jobs", "-j", min=1, help="claims run in parallel")] = 4,
    max_runs: MaxRunsOpt = None,
) -> None:
    """Reproduce the benchmark: Kang 2018 real-data claims and planted errors on PBMC3k.

    Downloads the datasets on first use and rewrites the committed bench/results.json and
    bench/results_table.md.
    """
    from recluster.bench.runner import run_bench
    from recluster.config import GridConfig

    parts = ("kang", "planted") if part is BenchPart.all else (part.value,)
    grid = GridConfig(max_runs=max_runs) if max_runs is not None else None
    run_bench(data_dir=data_dir, out_dir=out, parts=parts, jobs=jobs, grid=grid)


@app.command()
def demo(
    data_dir: Annotated[Path, typer.Option("--data-dir", help="where Kang 2018 is cached")] = Path(
        "data"
    ),
    out: Annotated[Path, typer.Option("--out", help="output directory")] = Path("results/demo"),
    jobs: Annotated[
        int | None,
        typer.Option(
            "--jobs",
            "-j",
            min=1,
            help="worker processes (default: cores - 1, max 8)",
        ),
    ] = None,
    no_open: Annotated[bool, typer.Option("--no-open", help="do not open the reports")] = False,
) -> None:
    """Run two claims on a subset of Kang 2018 (downloaded on first use): one true, one false
    by construction; opens both reports.
    """
    from recluster.demo import run_demo

    results = run_demo(data_dir=data_dir, out_dir=out, jobs=jobs, open_reports=not no_open)
    if any(r.verdict != r.expected for r in results):
        typer.secho("a demo claim did not get its expected verdict", err=True, fg="red")
        raise typer.Exit(1)


@app.command()
def parse(
    claim: Annotated[str, typer.Argument(help="the claim in plain English, in quotes")],
    data: DataOpt,
    out: Annotated[Path, typer.Option("--out", help="file to write the claim YAML to")] = Path(
        "claim.yaml"
    ),
    model: Annotated[str | None, typer.Option("--model", help="Anthropic model id")] = None,
    yes: Annotated[bool, typer.Option("--yes", "-y", help="skip the confirmation prompt")] = False,
    run: Annotated[bool, typer.Option("--run", help="run the claim after confirmation")] = False,
    out_dir: Annotated[
        Path, typer.Option("--out-dir", help="results folder, used with --run")
    ] = Path("results"),
    audit: Annotated[
        bool,
        typer.Option("--audit", help="the data is a processed object (same as `recluster audit`)"),
    ] = False,
) -> None:
    """Turn a plain-English claim into a claim YAML with an LLM, checked against the data.

    Optional: needs ANTHROPIC_API_KEY. Without it, copy examples/claims/TEMPLATE.yaml and
    check it with `recluster check`.
    """
    from recluster.data import DataError, load_adata
    from recluster.llm.common import LLMUnavailable, get_client
    from recluster.llm.parse import parse_claim, yaml_with_header
    from recluster.spec import SpecError

    try:
        client = get_client()
    except LLMUnavailable as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    try:
        adata = load_adata(data)
        result = parse_claim(claim, adata, client=client, model=model, audit=audit)
    except (LLMUnavailable, SpecError, DataError) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    text = yaml_with_header(result, claim)
    typer.echo(
        f"Proposed claim YAML (LLM output from {result.model}, {result.attempts} attempt(s)):\n"
    )
    typer.echo(text)
    if not yes and not typer.confirm(f"Write this claim YAML to {out}?", default=False):
        typer.echo("not written; edit by hand or rephrase the claim")
        raise typer.Exit(1)
    out.write_text(text)
    typer.echo(f"wrote {out}")
    if run:
        _execute(data, out, out_dir, None, None, None, None, False, audit=audit)


@app.command()
def summarize(
    provenance: Annotated[
        Path,
        typer.Argument(
            help="provenance.json from a results folder, "
            "e.g. `results/<claim-name>/provenance.json`",
            exists=True,
        ),
    ],
    model: Annotated[str | None, typer.Option("--model", help="Anthropic model id")] = None,
) -> None:
    """Write a plain-English LLM summary of an existing run's computed results (optional)."""
    import json

    from recluster.llm.common import LLMUnavailable
    from recluster.llm.summarize import SummaryRejected, summarize_provenance, summary_markdown

    try:
        text, model_id = summarize_provenance(json.loads(provenance.read_text()), model=model)
    except (LLMUnavailable, SummaryRejected) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    target = provenance.parent / "llm_summary.md"
    target.write_text(summary_markdown(text, model_id))
    typer.echo(f"LLM summary of computed results ({model_id}):\n{text}\n\nwrote {target}")


@app.command()
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


if __name__ == "__main__":  # pragma: no cover
    app()
