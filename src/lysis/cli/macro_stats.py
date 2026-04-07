"""``lysis macro-stats`` — print macroscale summary statistics for a simulation Run."""

import os

import click

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Rich display metadata
# ---------------------------------------------------------------------------

# Abbreviated column headers for the Rich table (may contain \n for multi-line)
_METRIC_SHORT = {
    "Degradation rate (%/min)": "Deg. Rate\n(%/min)",
    "Lysis lag time (min)": "Lag Time\n(min)",
    "Time to full clot degradation (min)": "Full Deg.\nTime (min)",
    "Percent of molecules that reached the back row": "Mol. Reach\nBack (%)",
    "First passage time (min)": "FPT\n(min)",
    "Front Velocity (microns/min)": "Front Vel.\n(\u03bcm/min)",
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="macro-stats")
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
@click.option(
    "--no-progress",
    is_flag=True,
    default=False,
    help="Suppress progress indicators.",
)
@click.option(
    "--markdown",
    "markdown_out",
    type=str,
    default=None,
    metavar="FILE",
    help=(
        "Output results as Markdown.  Use '-' to print to the console, "
        "or provide a filename to write to a file.  "
        "Progress indicators are suppressed automatically.  "
        "Example: --markdown -, --markdown report.md"
    ),
)
@click.pass_context
def macro_stats(ctx, path, sort_mode, no_progress, markdown_out):
    """Print macroscale summary statistics for one or more simulation Runs.

    PATH may be a single HDF5 file or a directory. When a directory is given,
    statistics for every .h5 file found in that directory are computed and
    displayed as a table, with rows ordered by --sort.

    Single-file Markdown output combines Mean and Std Dev into a single
    ``mean ± std`` value column (unlike the previous three-column format).

    \b
    Examples:
        lysis macro-stats data/TB-xi__1_582_867.h5
        lysis macro-stats data/lysis-front/
        lysis macro-stats data/lysis-front/ --sort alpha
        lysis macro-stats data/ --markdown -
        lysis macro-stats data/ --markdown report.md
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import macro_stats_table
    from lysis.cli.display import emit_markdown, stats_df_to_markdown, stats_df_to_rich
    from lysis.tools.runcode_sort import smart_sort
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
    )

    console = ctx.obj["console"]
    path = os.path.abspath(path)

    # Suppress progress when producing structured markdown output
    if markdown_out is not None:
        no_progress = True

    # -----------------------------------------------------------------------
    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)

        if not no_progress:
            with console.status(f"Computing statistics for {run_code}..."):
                stats = _load_run_stats(data_root, run_code, console)
        else:
            stats = _load_run_stats(data_root, run_code, console)

        if stats is None:
            ctx.exit(1)
            return

        rows = {run_code: stats}
        df = macro_stats_table(rows)

        if markdown_out is not None:
            emit_markdown(
                stats_df_to_markdown(df, "Run", single_code=run_code),
                markdown_out,
                console,
            )
        else:
            console.print(f"[bold]{run_code}[/bold]\n")
            for col, val in df.loc[run_code].items():
                console.print(f"  {col}: {val}")

    # -----------------------------------------------------------------------
    else:
        # --- directory mode ---
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

        _pctx = (
            Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeRemainingColumn(),
                console=console,
            )
            if not no_progress
            else nullcontext()
        )

        with _pctx as prog:
            _task = (
                prog.add_task("Computing macro statistics...", total=len(run_codes))
                if prog is not None
                else None
            )
            for run_code in run_codes:
                stats = _load_run_stats(path, run_code, console)
                if prog is not None:
                    prog.advance(_task)
                if stats is not None:
                    rows[run_code] = stats

        if not rows:
            ctx.exit(1)
            return

        # Preserve sort order, skipping any failed runs
        ordered = {rc: rows[rc] for rc in run_codes if rc in rows}
        df = macro_stats_table(ordered)

        if markdown_out is not None:
            emit_markdown(
                stats_df_to_markdown(df, "Run"),
                markdown_out,
                console,
            )
        else:
            console.print(stats_df_to_rich(df, "Run", _METRIC_SHORT))
