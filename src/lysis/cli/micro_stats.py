"""``lysis micro-stats`` — print microscale statistics for a simulation Run.

Mirrors the "## Microscale Measures" notebook section (cell 1.9.1), which
reports fiber lysis counts and times together with tPA leaving times across
all microscale simulations.
"""

import os

import click

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Rich display metadata
# ---------------------------------------------------------------------------

# Abbreviated column headers for the Rich table (may contain \n for multi-line)
_METRIC_SHORT = {
    "Fibers Degraded": "Fibers\nDegraded",
    "Mean Lysis Time (min)": "Mean Lysis\nTime (min)",
    "Median Lysis Time (min)": "Median Lysis\nTime (min)",
    "Mean tPA Leaving Time (sec)": "Mean tPA\nLeaving (sec)",
    "Median tPA Leaving Time (sec)": "Median tPA\nLeaving (sec)",
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_run_micro_stats(data_root, run_code, console):
    """Load a Run and compute microscale statistics.

    Delegates to :func:`~lysis.analysis.microscale.compute_micro_statistics`.

    :param data_root: Directory containing the HDF5 file.
    :type data_root: str
    :param run_code: Run code (HDF5 filename without extension).
    :type run_code: str
    :param console: Rich Console for error messages.
    :return: Stats dict with keys ``fibers_degraded``, ``lysis_time_mean``,
        ``lysis_time_std``, ``lysis_time_median``, ``tpa_leaving_mean``,
        ``tpa_leaving_std``, ``tpa_leaving_median``; or ``None`` on error.
    :rtype: dict or None
    """
    from lysis.analysis.microscale import compute_micro_statistics
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None

    try:
        return compute_micro_statistics(run)
    except Exception as e:
        console.print(f"[red]Error computing statistics for {run_code}:[/red] {e}")
        return None
    finally:
        run.data.close()


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="micro-stats")
@click.argument("path", type=click.Path(exists=True, file_okay=True, dir_okay=True))
@click.option(
    "--sort",
    "sort_mode",
    type=click.Choice(["smart", "alpha"], case_sensitive=False),
    default="smart",
    show_default=True,
    help=(
        "Sort order for directory mode. "
        "'smart' detects roman numerals and python-style integers and sorts "
        "them numerically; 'alpha' uses plain lexicographic order."
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
        "Example: --markdown -, --markdown stats.md"
    ),
)
@click.pass_context
def micro_stats(ctx, path, sort_mode, no_progress, markdown_out):
    """Print microscale statistics for one or more simulation Runs.

    Reports fiber degradation counts and lysis times together with tPA leaving
    times across all microscale simulations.  Mirrors the notebook's
    "## Microscale Measures" section.

    PATH may be a single HDF5 file or a directory.  In directory mode every
    .h5 file is processed and results are shown as a table with one row per
    Run.

    \b
    Examples:
        lysis micro-stats data/TB-xi__1_582_867.h5
        lysis micro-stats data/lysis-front/
        lysis micro-stats data/lysis-front/ --sort alpha
        lysis micro-stats data/ --markdown -
        lysis micro-stats data/ --markdown stats.md
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import micro_stats_table
    from lysis.tools.display import emit_markdown, stats_df_to_markdown, stats_df_to_rich
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
            with console.status(f"Computing microscale statistics for {run_code}..."):
                stats = _load_run_micro_stats(data_root, run_code, console)
        else:
            stats = _load_run_micro_stats(data_root, run_code, console)

        if stats is None:
            ctx.exit(1)
            return

        rows = {run_code: stats}
        df = micro_stats_table(rows)

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
                prog.add_task("Computing micro statistics...", total=len(run_codes))
                if prog is not None
                else None
            )
            for run_code in run_codes:
                stats = _load_run_micro_stats(path, run_code, console)
                if prog is not None:
                    prog.advance(_task)
                if stats is not None:
                    rows[run_code] = stats

        if not rows:
            ctx.exit(1)
            return

        # Preserve sort order, skipping any failed runs
        ordered = {rc: rows[rc] for rc in run_codes if rc in rows}
        df = micro_stats_table(ordered)

        if markdown_out is not None:
            emit_markdown(
                stats_df_to_markdown(df, "Run"),
                markdown_out,
                console,
            )
        else:
            console.print(stats_df_to_rich(df, "Run", _METRIC_SHORT))
