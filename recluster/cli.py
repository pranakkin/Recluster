"""Command-line interface."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

import typer

from recluster import __version__

app = typer.Typer(
    name="recluster",
    help="An independent stress test for single-cell RNA-seq claims.",
    no_args_is_help=True,
    add_completion=False,
)

DataOpt = Annotated[
    Path, typer.Option("--data", help="AnnData .h5ad file", exists=True, dir_okay=False)
]
ClaimOpt = Annotated[
    Path, typer.Option("--claim", help="claim spec YAML", exists=True, dir_okay=False)
]
OutOpt = Annotated[Path, typer.Option("--out", help="output directory")]
GridOpt = Annotated[
    Path | None, typer.Option("--grid", help="grid config YAML (overrides defaults)")
]
ThrOpt = Annotated[
    Path | None, typer.Option("--thresholds", help="thresholds YAML (overrides defaults)")
]
MaxRunsOpt = Annotated[int | None, typer.Option("--max-runs", help="cap the grid (quick look)")]
RandomOpt = Annotated[int | None, typer.Option("--n-random", help="number of random combinations")]
QuietOpt = Annotated[bool, typer.Option("--quiet", "-q", help="no per-run progress")]
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
        "--summarize", help="add an LLM summary of the computed results (needs an API key)"
    ),
]


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
    from recluster.config import GridConfig, Thresholds
    from recluster.data import DataError
    from recluster.engine import run_claim
    from recluster.spec import ClaimSpec, SpecError

    try:
        spec = ClaimSpec.from_yaml(claim)
    except Exception as exc:  # pydantic errors are verbose; show them plainly
        typer.secho(f"invalid claim spec: {exc}", err=True, fg="red")
        raise typer.Exit(2) from None
    cfg = GridConfig.from_yaml(grid) if grid else GridConfig()
    if max_runs is not None:
        cfg = cfg.model_copy(update={"max_runs": max_runs})
    if n_random is not None:
        cfg = cfg.model_copy(update={"n_random": n_random})
    thr = Thresholds.from_yaml(thresholds) if thresholds else Thresholds()
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
    """Stress-test a claim against raw counts, re-deriving everything from scratch."""
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
    """Audit a processed object: re-derive from its stored raw counts and check its labels."""
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


@app.command()
def bench(
    data_dir: Annotated[
        Path, typer.Option("--data-dir", help="where datasets are downloaded")
    ] = Path("data"),
    out: Annotated[Path, typer.Option("--out", help="where the benchmark reports go")] = Path(
        "bench/out"
    ),
    part: Annotated[str, typer.Option("--part", help="kang | planted | all")] = "all",
    jobs: Annotated[int, typer.Option("--jobs", help="claims run in parallel")] = 4,
    max_runs: MaxRunsOpt = None,
) -> None:
    """Reproduce the benchmark: Kang 2018 real-data claims and planted errors on PBMC3k."""
    from recluster.bench.runner import run_bench
    from recluster.config import GridConfig

    if part not in ("kang", "planted", "all"):
        typer.secho("--part must be kang, planted or all", err=True, fg="red")
        raise typer.Exit(2)
    parts = ("kang", "planted") if part == "all" else (part,)
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
    """Two claims on a subset of Kang 2018: one true, one false by construction."""
    from recluster.demo import run_demo

    results = run_demo(data_dir=data_dir, out_dir=out, jobs=jobs, open_reports=not no_open)
    if any(r.verdict != r.expected for r in results):
        typer.secho("a demo claim did not get its expected verdict", err=True, fg="red")
        raise typer.Exit(1)


@app.command()
def parse(
    claim: Annotated[str, typer.Argument(help="the claim in plain English, in quotes")],
    data: DataOpt,
    out: Annotated[Path, typer.Option("--out", help="where to write the proposed YAML")] = Path(
        "claim.yaml"
    ),
    model: Annotated[str | None, typer.Option("--model", help="Anthropic model id")] = None,
    yes: Annotated[bool, typer.Option("--yes", "-y", help="skip the confirmation prompt")] = False,
    run: Annotated[bool, typer.Option("--run", help="run the claim after confirmation")] = False,
    out_dir: Annotated[Path, typer.Option("--out-dir", help="results directory when --run")] = Path(
        "results"
    ),
    audit: Annotated[
        bool, typer.Option("--audit", help="treat the data as a processed object")
    ] = False,
) -> None:
    """Turn a plain-English claim into a claim spec with an LLM, checked against the data.

    Optional: needs ANTHROPIC_API_KEY. Without it, write the YAML by hand (see README).
    """
    from recluster.data import load_adata
    from recluster.llm.common import LLMUnavailable, get_client
    from recluster.llm.parse import parse_claim, yaml_with_header
    from recluster.spec import SpecError

    try:
        client = get_client()
    except LLMUnavailable as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    adata = load_adata(data)
    try:
        result = parse_claim(claim, adata, client=client, model=model, audit=audit)
    except (LLMUnavailable, SpecError) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    text = yaml_with_header(result, claim)
    typer.echo(
        f"Proposed claim spec (LLM output from {result.model}, {result.attempts} attempt(s)):\n"
    )
    typer.echo(text)
    if not yes and not typer.confirm(f"Write this spec to {out}?", default=False):
        typer.echo("not written; edit by hand or rephrase the claim")
        raise typer.Exit(1)
    out.write_text(text)
    typer.echo(f"wrote {out}")
    if run:
        _execute(data, out, out_dir, None, None, None, None, False, audit=audit)


@app.command()
def summarize(
    provenance: Annotated[
        Path, typer.Argument(help="provenance.json of a finished run", exists=True)
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
