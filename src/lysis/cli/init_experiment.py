"""``lysis init-experiment`` — initialise an Experiment from a parameter CSV."""

import os
import sys

import click

from lysis.cli import cli


@cli.command("init-experiment")
@click.argument("csv_path", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("data_root", type=click.Path(file_okay=False, path_type=str))
@click.option(
    "--name",
    default=None,
    show_default=False,
    help=(
        "Experiment name (used as the folder name under DATA_ROOT). "
        "Defaults to the CSV filename stem."
    ),
)
@click.option(
    "--description",
    default="",
    show_default=False,
    help="Optional prose description stored in experiment.json.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help=(
        "Validate and resolve all parameters but do not create any files or folders. "
        "Prints a summary of what would be created."
    ),
)
@click.option(
    "--no-progress",
    is_flag=True,
    default=False,
    help="Suppress progress indicators.",
)
@click.pass_context
def init_experiment(ctx, csv_path, data_root, name, description, dry_run, no_progress):
    """Initialise an Experiment from a parameter CSV file.

    Reads CSV_PATH (one column = one Run), resolves any missing dependent
    parameters algebraically, validates consistency, then creates the
    experiment folder under DATA_ROOT with an HDF5 file for each Run.

    \b
    Examples:
        lysis init-experiment runs.csv /data/experiments/
        lysis init-experiment runs.csv /data/experiments/ --name fiber-sweep
        lysis init-experiment runs.csv /data/experiments/ --dry-run
    """
    from contextlib import nullcontext

    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
    )
    from rich.table import Table

    from lysis.config.experiment import Experiment
    from lysis.config.param_resolver import ParameterConflict, UnderdeterminedParameters

    console = ctx.obj["console"]

    # Ensure data_root exists (or create it silently)
    os.makedirs(data_root, exist_ok=True)

    # ── Dry-run path ──────────────────────────────────────────────────────
    if dry_run:
        try:
            if not no_progress:
                with console.status("Validating parameters..."):
                    exp = Experiment.from_csv(
                        csv_path, data_root, name=name, description=description, dry_run=True
                    )
            else:
                exp = Experiment.from_csv(
                    csv_path, data_root, name=name, description=description, dry_run=True
                )
        except (ParameterConflict, UnderdeterminedParameters, ValueError) as exc:
            _print_error(console, csv_path, exc)
            ctx.exit(1)
            return

        console.print(
            f"[bold green]Dry run OK[/bold green] — {len(exp.runs)} run(s) validated, "
            "no files written."
        )
        _print_run_table(console, exp, dry_run=True)
        return

    # ── Real run ──────────────────────────────────────────────────────────
    try:
        exp = Experiment.from_csv(
            csv_path, data_root, name=name, description=description, dry_run=False
        )
    except FileExistsError as exc:
        stem = os.path.splitext(os.path.basename(csv_path))[0]
        console.print(f"[bold red]Error:[/bold red] {exc}")
        console.print(
            f"Tip: use [bold]--name[/bold] to choose a different folder name "
            f"(e.g. [bold]--name {stem}-v2[/bold])."
        )
        ctx.exit(1)
        return
    except (ParameterConflict, UnderdeterminedParameters, ValueError) as exc:
        _print_error(console, csv_path, exc)
        ctx.exit(1)
        return

    console.print(
        f"[bold green]Experiment created:[/bold green] {exp.path}"
    )
    _print_run_table(console, exp, dry_run=False)


def _print_error(console, csv_path, exc):
    """Print a formatted error block for parameter resolution failures."""
    from rich.panel import Panel

    title = f"[bold red]Parameter error in {os.path.basename(csv_path)}[/bold red]"
    console.print(Panel(str(exc), title=title, border_style="red"))


def _print_run_table(console, exp, *, dry_run: bool):
    """Print a summary table of run codes and descriptions."""
    from rich.table import Table

    tbl = Table(show_header=True, header_style="bold")
    tbl.add_column("Run code")
    tbl.add_column("Description")
    if not dry_run:
        tbl.add_column("HDF5")

    for run in exp.runs:
        desc = getattr(run, "_description", "") or ""
        if dry_run:
            tbl.add_row(run.run_code, desc)
        else:
            h5_name = f"{run.run_code}.h5"
            tbl.add_row(run.run_code, desc, h5_name)

    console.print(tbl)
