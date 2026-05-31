"""``lysis init-experiment`` — initialise an Experiment from a parameter CSV."""

import os
import sys

import click

from lysis.cli import cli
from lysis.cli._provenance import allow_dirty_option, enforce_lysis_clean


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
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help=(
        "If the experiment folder already exists, delete it and recreate from scratch. "
        "WARNING: all existing data in the folder will be lost."
    ),
)
@click.option(
    "--random-entropy",
    is_flag=True,
    default=False,
    help=(
        "Draw a fresh random micro_seed/macro_seed for every Run, ignoring any "
        "value in the CSV. By default a fresh value is drawn only for blank seed "
        "cells; an explicit seed is always respected."
    ),
)
@allow_dirty_option
@click.pass_context
def init_experiment(
    ctx, csv_path, data_root, name, description, dry_run, no_progress, force,
    random_entropy, allow_dirty,
):
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

    # Gate: refuse to write provenance if src/lysis/ is dirty (unless overridden).
    enforce_lysis_clean(ctx, allow_dirty)

    # Ensure data_root exists (or create it silently)
    os.makedirs(data_root, exist_ok=True)

    # ── Dry-run path ──────────────────────────────────────────────────────
    if dry_run:
        try:
            if not no_progress:
                with console.status("Validating parameters..."):
                    exp = Experiment.from_csv(
                        csv_path, data_root, name=name, description=description,
                        dry_run=True, random_entropy=random_entropy,
                    )
            else:
                exp = Experiment.from_csv(
                    csv_path, data_root, name=name, description=description,
                    dry_run=True, random_entropy=random_entropy,
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
    # If --force, delete the existing experiment folder before recreating.
    if force:
        import shutil

        stem = os.path.splitext(os.path.basename(csv_path))[0]
        exp_name = name or stem
        exp_path = os.path.join(data_root, exp_name)
        if os.path.isdir(exp_path):
            shutil.rmtree(exp_path)
            console.print(
                f"[yellow]Removed existing folder:[/yellow] {exp_path}"
            )

    try:
        exp = Experiment.from_csv(
            csv_path, data_root, name=name, description=description,
            dry_run=False, random_entropy=random_entropy,
        )
    except FileExistsError as exc:
        stem = os.path.splitext(os.path.basename(csv_path))[0]
        console.print(f"[bold red]Error:[/bold red] {exc}")
        console.print(
            f"Tip: use [bold]--name[/bold] to choose a different folder name "
            f"(e.g. [bold]--name {stem}-v2[/bold]), or use [bold]--force[/bold] "
            f"to replace the existing folder."
        )
        ctx.exit(1)
        return
    except (ParameterConflict, UnderdeterminedParameters, ValueError) as exc:
        _print_error(console, csv_path, exc)
        ctx.exit(1)
        return

    # Microscale init provenance is stamped inside DataStore.create()
    # (reached via Experiment.from_csv), so no separate stamp is needed here.

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
