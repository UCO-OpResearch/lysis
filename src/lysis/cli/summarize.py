"""``lysis summarize`` — print summary statistics for a simulation Run."""

import os

import click

from lysis.cli import cli

_METRIC_SHORT = {
    "Degradation rate (%/min)": "Deg. Rate\n(%/min)",
    "Lysis lag time (min)": "Lag Time\n(min)",
    "Time to full clot degradation (min)": "Full Deg.\nTime (min)",
    "Percent of molecules that reached the back row": "Mol. Reach\nBack (%)",
    "First passage time (min)": "FPT\n(min)",
    "Front Velocity (microns/min)": "Front Vel.\n(\u03bcm/min)",
}


def _load_run_stats(data_root, run_code, console):
    """Load a Run and return its statistics, or None on error."""
    from lysis.analysis.degradation import compute_run_statistics
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.macro_params = run.data.macro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None

    try:
        return compute_run_statistics(run)
    except Exception as e:
        console.print(f"[red]Error computing statistics for {run_code}:[/red] {e}")
        return None
    finally:
        run.data.close()


@cli.command()
@click.argument("path", type=click.Path(exists=True, file_okay=True, dir_okay=True))
@click.option(
    "--sort",
    "sort_mode",
    type=click.Choice(["smart", "alpha"], case_sensitive=False),
    default="smart",
    show_default=True,
    help=(
        "Sort order for directory mode. "
        "'smart' detects roman numerals and python-style integers and sorts them "
        "numerically; 'alpha' uses plain lexicographic order."
    ),
)
@click.pass_context
def summarize(ctx, path, sort_mode):
    """Print summary statistics for one or more simulation Runs.

    PATH may be a single HDF5 file or a directory. When a directory is given,
    statistics for every .h5 file found in that directory are computed and
    displayed as a table, with rows ordered by --sort.

    \b
    Examples:
        lysis summarize data/TB-xi__1_582_867.h5
        lysis summarize data/lysis-front/
        lysis summarize data/lysis-front/ --sort alpha
    """
    console = ctx.obj["console"]
    path = os.path.abspath(path)

    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)
        stats = _load_run_stats(data_root, run_code, console)
        if stats is None:
            ctx.exit(1)
            return
        console.print(f"[bold]{run_code}[/bold]\n")
        for metric in stats.index.get_level_values(0).unique():
            x = stats[metric]
            console.print(
                f"  {metric}: {x['Mean']:,.3f} \u00b1 {x['Standard Deviation']:,.3f}"
            )
    else:
        # --- directory mode: table of all .h5 files ---
        from lysis.tools.runcode_sort import smart_sort
        from rich.table import Table

        h5_files = [f for f in os.listdir(path) if f.lower().endswith(".h5")]
        if not h5_files:
            console.print(f"[yellow]No .h5 files found in {path}[/yellow]")
            ctx.exit(1)
            return

        run_codes = [os.path.splitext(f)[0] for f in h5_files]
        if sort_mode == "smart":
            run_codes = smart_sort(run_codes)
        else:
            run_codes = sorted(run_codes)

        rows = {}
        first_stats = None
        for run_code in run_codes:
            stats = _load_run_stats(path, run_code, console)
            if stats is not None:
                rows[run_code] = stats
                if first_stats is None:
                    first_stats = stats

        if not rows:
            ctx.exit(1)
            return

        metrics = first_stats.index.get_level_values(0).unique().tolist()

        table = Table(show_header=True, header_style="bold", show_lines=True)
        table.add_column("Run", style="bold", no_wrap=True)
        for metric in metrics:
            table.add_column(_METRIC_SHORT.get(metric, metric), justify="right")

        for run_code, stats in rows.items():
            cells = [run_code]
            for metric in metrics:
                x = stats[metric]
                cells.append(
                    f"{x['Mean']:,.3f}\n\u00b1 {x['Standard Deviation']:,.3f}"
                )
            table.add_row(*cells)

        console.print(table)
