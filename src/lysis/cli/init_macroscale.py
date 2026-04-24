"""``lysis init-macroscale`` — initialise macroscale structure in an HDF5 file.

Accepts either:

* A **folder path** (experiment folder containing ``experiment.json``) — reads
  macro parameters from the JSON and initialises every run in the experiment.
* An **HDF5 file path** (``*.h5``) — initialises that single run.  Macro
  parameters are constructed from defaults; use ``--param NAME=VALUE`` to
  override individual values.

In both cases the microscale datasets must already contain completed simulation
results so that ``forced_unbind`` can be computed automatically.
"""

import os
import sys

import click

from lysis.cli import cli


def _coerce_value(val_str: str):
    """Convert a string to int, float, bool, or leave as-is for Pint parsing."""
    stripped = val_str.strip()
    if stripped.lower() == "true":
        return True
    if stripped.lower() == "false":
        return False
    try:
        return int(stripped)
    except ValueError:
        pass
    try:
        return float(stripped)
    except ValueError:
        pass
    return stripped


@cli.command("init-macroscale")
@click.argument("path", type=click.Path(exists=True, path_type=str))
@click.option(
    "--param",
    "params",
    multiple=True,
    metavar="NAME=VALUE",
    help=(
        "Override a MacroParameters field (repeatable). "
        "Only used when PATH is an HDF5 file. "
        "Example: --param rows=121 --param cols=93"
    ),
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help=(
        "Validate that microscale data is present and compute forced_unbind, "
        "but do not write any macroscale structure to disk."
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
        "If macroscale structure already exists, replace it with blank tables. "
        "WARNING: any existing macroscale simulation data will be lost."
    ),
)
@click.pass_context
def init_macroscale(ctx, path, params, dry_run, no_progress, force):
    """Initialise macroscale structure in one or more HDF5 files.

    PATH may be an experiment folder (initialises all runs) or a single HDF5
    file (initialises that run only).  Microscale simulations must have
    completed before running this command; ``forced_unbind`` is computed
    automatically from the completed microscale output.

    \b
    Examples:
        lysis init-macroscale /data/experiments/fiber-sweep/
        lysis init-macroscale /data/experiments/fiber-sweep/run-00.h5
        lysis init-macroscale run.h5 --param rows=64 --param cols=64
        lysis init-macroscale /data/experiments/fiber-sweep/ --dry-run
    """
    from pathlib import Path

    from lysis.config.experiment import Experiment
    from lysis.config.paramcheck import load_macro_params
    from lysis.dataio.datastore import DataStore

    console = ctx.obj["console"]
    path_obj = Path(path)

    # ── Parse --param overrides ───────────────────────────────────────────
    param_overrides = {}
    for p in params:
        if "=" not in p:
            console.print(
                f"[bold red]Error:[/bold red] --param must be NAME=VALUE, "
                f"got: {p!r}"
            )
            ctx.exit(1)
            return
        name, _, value_str = p.partition("=")
        param_overrides[name.strip()] = _coerce_value(value_str.strip())

    # ── Dispatch: folder or single H5 file ───────────────────────────────
    if path_obj.is_dir():
        _init_experiment_folder(
            ctx, console, path_obj, param_overrides, dry_run, no_progress, force
        )
    elif path_obj.suffix == ".h5":
        _init_single_h5(
            ctx, console, path_obj, param_overrides, dry_run, no_progress, force
        )
    else:
        console.print(
            f"[bold red]Error:[/bold red] PATH must be an experiment folder "
            f"or an HDF5 file (.h5), got: {path!r}"
        )
        ctx.exit(1)


# ─── Folder mode ──────────────────────────────────────────────────────────────


def _init_experiment_folder(ctx, console, folder_path, param_overrides, dry_run, no_progress, force):
    """Initialise macroscale for every run in an experiment folder."""
    from lysis.config.experiment import Experiment
    from lysis.dataio.datastore import DataStore

    if param_overrides:
        console.print(
            "[yellow]Warning:[/yellow] --param overrides are ignored in folder "
            "mode; macro parameters are read from experiment.json."
        )

    # Load the Experiment (reads experiment.json + HDF5 micro_params)
    try:
        exp = Experiment.load(folder_path)
    except FileNotFoundError as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        ctx.exit(1)
        return
    except Exception as exc:
        console.print(f"[bold red]Error loading experiment:[/bold red] {exc}")
        ctx.exit(1)
        return

    errors = []
    results = []

    for run in exp.runs:
        h5_path = os.path.join(exp.path, f"{run.run_code}.h5")
        try:
            if dry_run:
                # Open read-only; just verify microscale data exists
                with DataStore(run.run_code, exp.path, mode="r") as ds:
                    _check_microscale_ready(ds, run.run_code)
                results.append((run.run_code, "dry-run OK"))
            else:
                with DataStore(run.run_code, exp.path, mode="a") as ds:
                    ds.initialize_macroscale(run.macro_params, force=force)
                results.append((run.run_code, "initialized"))
        except Exception as exc:
            errors.append((run.run_code, str(exc)))

    _print_results(console, results, errors, dry_run)

    if errors:
        ctx.exit(1)


# ─── Single H5 mode ───────────────────────────────────────────────────────────


def _init_single_h5(ctx, console, h5_path, param_overrides, dry_run, no_progress, force):
    """Initialise macroscale for a single HDF5 file."""
    from lysis.config.paramcheck import load_macro_params
    from lysis.dataio.datastore import DataStore

    run_code = h5_path.stem
    run_dir = str(h5_path.parent)

    # Read micro_params from the file (required to construct MacroParameters)
    try:
        with DataStore(run_code, run_dir, mode="r") as ds:
            micro_params = ds.micro_params
            if micro_params is None:
                console.print(
                    f"[bold red]Error:[/bold red] No microscale parameters found "
                    f"in {h5_path.name}. Is this a valid lysis HDF5 file?"
                )
                ctx.exit(1)
                return
            if dry_run:
                _check_microscale_ready(ds, run_code)
    except Exception as exc:
        console.print(f"[bold red]Error reading {h5_path.name}:[/bold red] {exc}")
        ctx.exit(1)
        return

    # Construct MacroParameters from defaults + overrides
    try:
        macro_params = load_macro_params({}, micro_params, overrides=param_overrides)
    except ValueError as exc:
        console.print(f"[bold red]Error in --param overrides:[/bold red] {exc}")
        ctx.exit(1)
        return

    if dry_run:
        console.print(
            f"[bold green]Dry run OK[/bold green] — "
            f"{run_code}: microscale data present, macroscale ready to initialize."
        )
        return

    # Initialize macroscale
    try:
        with DataStore(run_code, run_dir, mode="a") as ds:
            ds.initialize_macroscale(macro_params, force=force)
    except Exception as exc:
        console.print(f"[bold red]Error initializing {h5_path.name}:[/bold red] {exc}")
        ctx.exit(1)
        return

    console.print(
        f"[bold green]Initialized:[/bold green] {h5_path.name} — "
        f"macroscale structure written."
    )


# ─── Helpers ──────────────────────────────────────────────────────────────────


def _check_microscale_ready(ds, run_code):
    """Raise ValueError if the DataStore has no microscale data."""
    from lysis.dataio.dataspec import dataspec
    from lysis.dataio.datastore import COMPATIBLE_DATASPEC_VERSION

    if "microscale_out" not in ds.collections:
        raise ValueError(
            f"{run_code}: microscale_out collection is absent. "
            "Ensure the HDF5 file was created with DataStore.create()."
        )
    micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
    pli_loc = micro_spec.data["tpa_unbound_by_pli"].data_location
    import h5py
    with h5py.File(ds._hdf5_path, "r") as f:
        if f[pli_loc].shape[0] == 0:
            raise ValueError(
                f"{run_code}: microscale datasets are empty. "
                "Ensure all microscale Simulations have completed."
            )


def _print_results(console, results, errors, dry_run):
    """Print a summary table of per-run outcomes."""
    from rich.table import Table

    tbl = Table(show_header=True, header_style="bold")
    tbl.add_column("Run code")
    tbl.add_column("Status")

    for run_code, status in results:
        tbl.add_row(run_code, f"[green]{status}[/green]")
    for run_code, msg in errors:
        tbl.add_row(run_code, f"[red]ERROR: {msg}[/red]")

    console.print(tbl)

    if dry_run and not errors:
        console.print(
            f"[bold green]Dry run OK[/bold green] — "
            f"{len(results)} run(s) validated, no files written."
        )
    elif not errors:
        console.print(
            f"[bold green]Done:[/bold green] "
            f"{len(results)} run(s) initialized."
        )
    else:
        console.print(
            f"[bold red]{len(errors)} error(s)[/bold red] — "
            f"{len(results)} succeeded, {len(errors)} failed."
        )
