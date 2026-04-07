"""``lysis macro-stats`` — print macroscale summary statistics for a simulation Run."""

import os

import click

from lysis.cli import cli

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


def _summarize_to_markdown(rows_dict, metrics, single_file_code=None):
    """Build a Markdown string from summarize data.

    In single-file mode (*single_file_code* given) produces a heading and a
    three-column table (Metric / Mean / Std Dev).  In directory mode produces
    a flat table with one row per Run and one column per metric (mean ± std).

    :param rows_dict: Maps run code → stats Series.
    :type rows_dict: dict[str, pandas.Series]
    :param metrics: Ordered list of metric names.
    :type metrics: list[str]
    :param single_file_code: Run code when operating on a single file, or
        ``None`` for directory mode.
    :type single_file_code: str or None
    :return: Markdown text.
    :rtype: str
    """
    if single_file_code is not None:
        stats = rows_dict[single_file_code]
        table_rows = [
            [
                m,
                f"{stats[m]['Mean']:,.3f}",
                f"{stats[m]['Standard Deviation']:,.3f}",
            ]
            for m in metrics
        ]
        return "\n".join(
            [f"## {single_file_code}", "", _md_table(["Metric", "Mean", "Std Dev"], table_rows)]
        )
    else:
        headers = ["Run"] + metrics
        table_rows = [
            [rc]
            + [
                f"{rows_dict[rc][m]['Mean']:,.3f} \u00b1 {rows_dict[rc][m]['Standard Deviation']:,.3f}"
                for m in metrics
            ]
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

    \b
    Examples:
        lysis macro-stats data/TB-xi__1_582_867.h5
        lysis macro-stats data/lysis-front/
        lysis macro-stats data/lysis-front/ --sort alpha
        lysis macro-stats data/ --markdown -
        lysis macro-stats data/ --markdown report.md
    """
    console = ctx.obj["console"]
    path = os.path.abspath(path)

    # Suppress progress when producing structured markdown output
    if markdown_out is not None:
        no_progress = True

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

        metrics = stats.index.get_level_values(0).unique().tolist()

        if markdown_out is not None:
            _emit_markdown(
                _summarize_to_markdown({run_code: stats}, metrics, run_code),
                markdown_out,
                console,
            )
        else:
            console.print(f"[bold]{run_code}[/bold]\n")
            for metric in metrics:
                x = stats[metric]
                console.print(
                    f"  {metric}: {x['Mean']:,.3f} \u00b1 {x['Standard Deviation']:,.3f}"
                )

    else:
        # --- directory mode: table of all .h5 files ---
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
                    if first_stats is None:
                        first_stats = stats

        if not rows:
            ctx.exit(1)
            return

        metrics = first_stats.index.get_level_values(0).unique().tolist()

        if markdown_out is not None:
            _emit_markdown(
                _summarize_to_markdown(rows, metrics),
                markdown_out,
                console,
            )
        else:
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
