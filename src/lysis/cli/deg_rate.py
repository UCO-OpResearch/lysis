"""``lysis deg-rate`` — print degradation-rate tables for simulation Runs.

Mirrors the "## Degradation Rate Table" notebook cell, which computes the
mean ± standard deviation of degradation rates across macroscale simulations
for configurable degradation intervals (e.g. 20 %→80 %).
"""

import os

import click

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Default intervals (start%, end%) matching the notebook cell
_DEFAULT_INTERVALS = [(20, 80), (20, 50), (50, 80)]

# Display format for rates (shown in %/min)
_RATE_FMT = "{:.4f}"

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


def _deg_rate_to_markdown(rows_dict, intervals, single_file_code=None):
    """Build a Markdown string from degradation-rate data.

    In single-file mode (*single_file_code* given) produces a heading and a
    three-column table (Interval / Mean (%/min) / Std Dev (%/min)).  In
    directory mode produces a flat table with one row per Run.

    :param rows_dict: Maps run code → {interval_key: (mean, std)}.
    :type rows_dict: dict[str, dict[tuple, tuple]]
    :param intervals: Ordered list of (start, end) integer percent pairs.
    :type intervals: list[tuple[int, int]]
    :param single_file_code: Run code for single-file mode, or ``None``.
    :type single_file_code: str or None
    :return: Markdown text.
    :rtype: str
    """
    interval_labels = [f"{s}% to {e}% (%/min)" for s, e in intervals]

    if single_file_code is not None:
        stats = rows_dict[single_file_code]
        table_rows = [
            [
                f"{s}% to {e}%",
                _RATE_FMT.format(stats[(s, e)][0]),
                _RATE_FMT.format(stats[(s, e)][1]),
            ]
            for s, e in intervals
        ]
        return "\n".join(
            [
                f"## {single_file_code}",
                "",
                _md_table(["Interval", "Mean (%/min)", "Std Dev (%/min)"], table_rows),
            ]
        )
    else:
        headers = ["Run"] + interval_labels
        table_rows = [
            [rc]
            + [
                f"{rows_dict[rc][(s, e)][0]:.4f} \u00b1 {rows_dict[rc][(s, e)][1]:.4f}"
                for s, e in intervals
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
# Interval helpers
# ---------------------------------------------------------------------------


def _parse_interval(s, param_name="interval"):
    """Parse a ``'START-END'`` string into an ``(int, int)`` tuple.

    :param s: Interval string, e.g. ``'20-80'``.
    :type s: str
    :param param_name: Name used in error messages.
    :type param_name: str
    :return: ``(start, end)`` tuple of integer percentages.
    :rtype: tuple[int, int]
    :raises click.BadParameter: If the format or range is invalid.
    """
    parts = s.split("-")
    if len(parts) != 2:
        raise click.BadParameter(
            f"Expected START-END (e.g. '20-80'), got: {s!r}", param_hint=param_name
        )
    try:
        start, end = int(parts[0]), int(parts[1])
    except ValueError:
        raise click.BadParameter(
            f"Interval endpoints must be integers, got: {s!r}", param_hint=param_name
        )
    if not (0 <= start < end <= 100):
        raise click.BadParameter(
            f"Interval must satisfy 0 ≤ start < end ≤ 100, got: {s!r}",
            param_hint=param_name,
        )
    return (start, end)


def _build_percent_markers(intervals):
    """Return a sorted list of fractional markers covering all interval endpoints.

    Always includes 0.0 and 1.0 as bookends.

    :param intervals: List of (start, end) integer percent pairs.
    :type intervals: list[tuple[int, int]]
    :return: Sorted list of fractions.
    :rtype: list[float]
    """
    endpoints = {0, 100}
    for start, end in intervals:
        endpoints.add(start)
        endpoints.add(end)
    return [p / 100 for p in sorted(endpoints)]


def _interval_label(start, end, with_units=True):
    """Human-readable label for an interval.

    :param start: Start percentage (integer).
    :param end: End percentage (integer).
    :param with_units: If True, append ``'\\n(%/min)'`` for Rich multi-line headers.
    :return: Label string.
    :rtype: str
    """
    label = f"{start}% to {end}%"
    return f"{label}\n(%/min)" if with_units else label


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_run_deg_rates(data_root, run_code, intervals, console):
    """Load a Run and compute degradation rates for each interval.

    :param data_root: Directory containing the HDF5 file.
    :type data_root: str
    :param run_code: Run code (HDF5 filename without extension).
    :type run_code: str
    :param intervals: List of (start, end) integer percent pairs.
    :type intervals: list[tuple[int, int]]
    :param console: Rich Console for error messages.
    :return: Dict mapping interval key → (mean_pct_per_min, std_pct_per_min),
        or ``None`` on error.
    :rtype: dict[tuple, tuple] or None
    """
    from lysis.analysis.degradation import degradation_rates
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.macro_params = run.data.macro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None

    try:
        slope_pairs = [(s / 100, e / 100) for s, e in intervals]
        pct_markers = _build_percent_markers(intervals)
        rates = degradation_rates(run, slope_pairs, pct_markers)
        # rates: shape (n_sims, n_intervals), fraction/min → multiply by 100 for %/min
        rates_pct = rates * 100
        return {
            ivl: (float(rates_pct[:, k].mean()), float(rates_pct[:, k].std()))
            for k, ivl in enumerate(intervals)
        }
    except Exception as e:
        console.print(f"[red]Error computing rates for {run_code}:[/red] {e}")
        return None
    finally:
        run.data.close()


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="deg-rate")
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
        "Example: --markdown -, --markdown rates.md"
    ),
)
@click.option(
    "--add",
    "add_intervals",
    multiple=True,
    metavar="START-END",
    help=(
        "Add a degradation interval to the table (repeatable).  "
        "START and END are integer percentages.  "
        "Example: --add 0-100 --add 0-50"
    ),
)
@click.option(
    "--drop",
    "drop_intervals",
    multiple=True,
    metavar="START-END",
    help=(
        "Remove a default interval from the table (repeatable).  "
        "Example: --drop 20-80"
    ),
)
@click.pass_context
def deg_rate(ctx, path, sort_mode, no_progress, markdown_out, add_intervals, drop_intervals):
    """Print degradation-rate tables for one or more simulation Runs.

    Computes the mean and standard deviation of the degradation rate (%%/min)
    across all macroscale simulations for each requested interval.  Mirrors
    the notebook's "## Degradation Rate Table".

    PATH may be a single HDF5 file or a directory.  In directory mode every
    .h5 file is processed and results are shown as a table with one row per
    Run and one column per interval.

    Default intervals: 20%%→80%%, 20%%→50%%, 50%%→80%%.

    \b
    Examples:
        lysis deg-rate data/TB-xi__1_582_867.h5
        lysis deg-rate data/lysis-front/
        lysis deg-rate data/ --add 0-100 --drop 20-80
        lysis deg-rate data/ --markdown -
        lysis deg-rate data/ --markdown rates.md
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

    # Build effective interval list
    drop_set = set()
    for raw in drop_intervals:
        try:
            drop_set.add(_parse_interval(raw))
        except click.BadParameter as e:
            console.print(f"[red]Error:[/red] --drop: {e.format_message()}")
            ctx.exit(1)
            return

    intervals = [ivl for ivl in _DEFAULT_INTERVALS if ivl not in drop_set]

    for raw in add_intervals:
        try:
            ivl = _parse_interval(raw)
        except click.BadParameter as e:
            console.print(f"[red]Error:[/red] --add: {e.format_message()}")
            ctx.exit(1)
            return
        if ivl not in intervals:
            intervals.append(ivl)

    if not intervals:
        console.print("[red]Error:[/red] No intervals remain after applying --drop.")
        ctx.exit(1)
        return

    # -----------------------------------------------------------------------
    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)

        if not no_progress:
            with console.status(f"Computing degradation rates for {run_code}..."):
                stats = _load_run_deg_rates(data_root, run_code, intervals, console)
        else:
            stats = _load_run_deg_rates(data_root, run_code, intervals, console)

        if stats is None:
            ctx.exit(1)
            return

        if markdown_out is not None:
            _emit_markdown(
                _deg_rate_to_markdown({run_code: stats}, intervals, run_code),
                markdown_out,
                console,
            )
        else:
            console.print(f"[bold]{run_code}[/bold]\n")
            for s, e in intervals:
                mean, std = stats[(s, e)]
                console.print(
                    f"  {s}% to {e}%: "
                    f"{_RATE_FMT.format(mean)} \u00b1 {_RATE_FMT.format(std)} %/min"
                )

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
                prog.add_task("Computing degradation rates...", total=len(run_codes))
                if prog is not None
                else None
            )
            for run_code in run_codes:
                stats = _load_run_deg_rates(path, run_code, intervals, console)
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
                _deg_rate_to_markdown(rows, intervals),
                markdown_out,
                console,
            )
        else:
            table = Table(show_header=True, header_style="bold", show_lines=True)
            table.add_column("Run", style="bold", no_wrap=True)
            for s, e in intervals:
                table.add_column(_interval_label(s, e), justify="right")

            for run_code in ordered:
                stats = rows[run_code]
                cells = [run_code]
                for s, e in intervals:
                    mean, std = stats[(s, e)]
                    cells.append(
                        f"{_RATE_FMT.format(mean)}\n\u00b1 {_RATE_FMT.format(std)}"
                    )
                table.add_row(*cells)

            console.print(table)
