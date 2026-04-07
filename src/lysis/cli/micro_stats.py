"""``lysis micro-stats`` — print microscale statistics for a simulation Run.

Mirrors the "## Microscale Measures" notebook section (cell 1.9.1), which
reports fiber lysis counts and times together with tPA leaving times across
all microscale simulations.
"""

import os

import click

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Metric definitions
# ---------------------------------------------------------------------------

# Ordered display names and their short Rich-table headers
_METRICS = [
    "Fibers Degraded",
    "Mean Lysis Time (min)",
    "Median Lysis Time (min)",
    "Mean tPA Leaving Time (sec)",
    "Median tPA Leaving Time (sec)",
]

_METRIC_SHORT = {
    "Fibers Degraded": "Fibers\nDegraded",
    "Mean Lysis Time (min)": "Mean Lysis\nTime (min)",
    "Median Lysis Time (min)": "Median Lysis\nTime (min)",
    "Mean tPA Leaving Time (sec)": "Mean tPA\nLeaving (sec)",
    "Median tPA Leaving Time (sec)": "Median tPA\nLeaving (sec)",
}

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _fmt_metric(key, stats):
    """Format a single metric value from *stats* as a display string.

    :param key: Metric name from :data:`_METRICS`.
    :type key: str
    :param stats: Stats dict as returned by :func:`_load_run_micro_stats`.
    :type stats: dict
    :return: Formatted string.
    :rtype: str
    """
    if key == "Fibers Degraded":
        return f"{stats['fibers_degraded']:,}"
    if key == "Mean Lysis Time (min)":
        return f"{stats['lysis_time_mean']:.3f} \u00b1 {stats['lysis_time_std']:.3f}"
    if key == "Median Lysis Time (min)":
        return f"{stats['lysis_time_median']:.3f}"
    if key == "Mean tPA Leaving Time (sec)":
        return f"{stats['tpa_leaving_mean']:.3f} \u00b1 {stats['tpa_leaving_std']:.3f}"
    if key == "Median tPA Leaving Time (sec)":
        return f"{stats['tpa_leaving_median']:.3f}"
    raise KeyError(key)


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------


def _md_table(headers, rows):
    """Build a GitHub-flavored Markdown table string.

    :param headers: Column header strings.
    :type headers: list[str]
    :param rows: Data rows, each a list of cell strings.
    :type rows: list[list[str]]
    :return: GFM table as a string.
    :rtype: str
    """
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _micro_stats_to_markdown(rows_dict, single_file_code=None):
    """Build a Markdown string from microscale-stats data.

    In single-file mode (*single_file_code* given) produces a heading and a
    two-column table (Metric / Value).  In directory mode produces a flat
    table with one row per Run and one column per metric.

    :param rows_dict: Maps run code → stats dict.
    :type rows_dict: dict[str, dict]
    :param single_file_code: Run code for single-file mode, or ``None``.
    :type single_file_code: str or None
    :return: Markdown text.
    :rtype: str
    """
    if single_file_code is not None:
        stats = rows_dict[single_file_code]
        table_rows = [[m, _fmt_metric(m, stats)] for m in _METRICS]
        return "\n".join(
            [
                f"## {single_file_code}",
                "",
                _md_table(["Metric", "Value"], table_rows),
            ]
        )
    else:
        headers = ["Run"] + _METRICS
        table_rows = [
            [rc] + [_fmt_metric(m, rows_dict[rc]) for m in _METRICS]
            for rc in rows_dict
        ]
        return _md_table(headers, table_rows)


def _emit_markdown(md_text, markdown_out, console):
    """Write *md_text* to stdout or a file.

    :param md_text: Markdown content to write.
    :type md_text: str
    :param markdown_out: ``"-"`` for stdout, or a file path.
    :type markdown_out: str
    :param console: Rich Console for status messages.
    """
    if markdown_out == "-":
        print(md_text)
    else:
        try:
            with open(markdown_out, "w", encoding="utf-8") as fh:
                fh.write(md_text)
                if not md_text.endswith("\n"):
                    fh.write("\n")
            console.print(f"[green]Markdown written to {markdown_out}[/green]")
        except OSError as e:
            console.print(f"[red]Error:[/red] Could not write to {markdown_out}: {e}")


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

    from lysis.tools.runcode_sort import smart_sort
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
    )
    from rich.table import Table

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

        if markdown_out is not None:
            _emit_markdown(
                _micro_stats_to_markdown({run_code: stats}, run_code),
                markdown_out,
                console,
            )
        else:
            console.print(f"[bold]{run_code}[/bold]\n")
            for m in _METRICS:
                console.print(f"  {m}: {_fmt_metric(m, stats)}")

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

        ordered = [rc for rc in run_codes if rc in rows]

        if markdown_out is not None:
            _emit_markdown(
                _micro_stats_to_markdown(rows),
                markdown_out,
                console,
            )
        else:
            table = Table(show_header=True, header_style="bold", show_lines=True)
            table.add_column("Run", style="bold", no_wrap=True)
            for m in _METRICS:
                table.add_column(_METRIC_SHORT.get(m, m), justify="right")

            for run_code in ordered:
                cells = [run_code] + [_fmt_metric(m, rows[run_code]) for m in _METRICS]
                table.add_row(*cells)

            console.print(table)
