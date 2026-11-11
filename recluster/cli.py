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
            progress=None if quiet else (lambda m: print(m, file=sys.stderr, flush=True)),
        )
    except (SpecError, DataError) as exc:
        typer.secho(str(exc), err=True, fg="red")
        raise typer.Exit(2) from None
    v = result.verdict
    typer.echo(f"{v.label}\t{v.n_holds}/{v.n_runs} runs\t{spec.name}")


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
) -> None:
    """Stress-test a claim against raw counts, re-deriving everything from scratch."""
    _execute(data, claim, out, grid, thresholds, max_runs, n_random, quiet, audit=False)


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
) -> None:
    """Audit a processed object: re-derive from its stored raw counts and check its labels."""
    _execute(data, claim, out, grid, thresholds, max_runs, n_random, quiet, audit=True)


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
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


if __name__ == "__main__":  # pragma: no cover
    app()
