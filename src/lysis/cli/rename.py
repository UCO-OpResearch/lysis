"""``lysis rename`` — rename an Experiment folder or a Run HDF5 file."""

import os

import click

from lysis.cli import cli


@cli.command("rename")
@click.argument("path", type=click.Path(exists=True, path_type=str))
@click.argument("new_name", type=str)
@click.pass_context
def rename(ctx, path, new_name):
    """Rename an Experiment folder or a Run HDF5 file.

    If PATH is a folder, the experiment folder is renamed and the
    ``name`` field in ``experiment.json`` is updated.

    If PATH is an HDF5 file (``.h5``), the file is renamed to
    ``NEW_NAME.h5`` and the previous run_code is recorded in the
    ``renamed_from`` root attribute.  If ``experiment.json`` is present
    in the same directory, the matching entry there is updated as well.

    NEW_NAME must be a simple identifier (no path separators).  For HDF5
    files a trailing ``.h5`` is stripped if supplied.

    \b
    Examples:
        lysis rename /data/experiments/old-name new-name
        lysis rename /data/experiments/exp/run01.h5 run01-rerun
    """
    console = ctx.obj["console"]

    if os.sep in new_name or (os.altsep and os.altsep in new_name):
        console.print(
            f"[red]Error:[/red] NEW_NAME must not contain path separators: "
            f"{new_name!r}"
        )
        ctx.exit(1)
        return

    if os.path.isdir(path):
        _rename_experiment(ctx, console, path, new_name)
    elif os.path.isfile(path) and path.endswith(".h5"):
        if new_name.endswith(".h5"):
            new_name = new_name[:-3]
        _rename_run(ctx, console, path, new_name)
    else:
        console.print(
            f"[red]Error:[/red] PATH must be an experiment folder or an "
            f"HDF5 (.h5) file: {path!r}"
        )
        ctx.exit(1)


def _rename_experiment(ctx, console, folder, new_name):
    from lysis.config.experiment import Experiment

    folder = os.path.abspath(folder)
    try:
        exp = Experiment.load(folder)
    except FileNotFoundError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        ctx.exit(1)
        return

    old_name = exp.name
    try:
        exp.rename(new_name)
    except (FileExistsError, ValueError) as exc:
        console.print(f"[red]Error:[/red] {exc}")
        ctx.exit(1)
        return

    console.print(
        f"[bold green]Renamed experiment[/bold green] "
        f"{old_name!r} → {new_name!r}\n  {exp.path}"
    )


def _rename_run(ctx, console, h5_path, new_run_code):
    from lysis.config.experiment import Experiment
    from lysis.config.run import Run

    h5_path = os.path.abspath(h5_path)
    folder = os.path.dirname(h5_path)
    old_run_code = os.path.splitext(os.path.basename(h5_path))[0]

    exp_json = os.path.join(folder, "experiment.json")
    if os.path.isfile(exp_json):
        try:
            exp = Experiment.load(folder)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            console.print(
                f"[red]Error:[/red] Failed to load experiment.json: {exc}"
            )
            ctx.exit(1)
            return
        try:
            exp.rename_run(old_run_code, new_run_code)
        except (KeyError, ValueError, FileExistsError, FileNotFoundError) as exc:
            console.print(f"[red]Error:[/red] {exc}")
            ctx.exit(1)
            return
        console.print(
            f"[bold green]Renamed run[/bold green] "
            f"{old_run_code!r} → {new_run_code!r} "
            f"(experiment.json updated)"
        )
        return

    run = Run(folder, old_run_code)
    try:
        run.rename(new_run_code)
    except (ValueError, FileExistsError, FileNotFoundError) as exc:
        console.print(f"[red]Error:[/red] {exc}")
        ctx.exit(1)
        return
    console.print(
        f"[bold green]Renamed run[/bold green] "
        f"{old_run_code!r} → {new_run_code!r}"
    )
