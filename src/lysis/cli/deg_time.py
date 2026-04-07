"""``lysis deg-time`` — print degradation-time tables for simulation Runs.

Mirrors the "## Degradation Time Table" notebook cell, which computes the
mean ± standard deviation of the time (minutes) for the clot to reach each
degradation milestone across macroscale simulations.
"""

import os

import click

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Default percentage milestones matching the notebook's degrade_percent_markers
_DEFAULT_MARKERS = [5, 20, 50, 80, 100]

# Display format for times (minutes)
_TIME_FMT = "{:.2f}"

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


def _deg_time_to_markdown(rows_dict, markers, single_file_code=None):
    """Build a Markdown string from degradation-time data.

    In single-file mode (*single_file_code* given) produces a heading and a
    three-column table (Milestone / Mean (min) / Std Dev (min)).  In
    directory mode produces a flat table with one row per Run.

    :param rows_dict: Maps run code → {marker_pct: (mean_min, std_min)}.
    :type rows_dict: dict[str, dict[int, tuple]]
    :param markers: Ordered list of integer percentage milestones.
    :type markers: list[int]
    :param single_file_code: Run code for single-file mode, or ``None``.
    :type single_file_code: str or None
    :return: Markdown text.
    :rtype: str
    """
    marker_labels = [f"{m}% (min)" for m in markers]

    if single_file_code is not None:
        stats = rows_dict[single_file_code]
        table_rows = [
            [
                f"{m}%",
                _TIME_FMT.format(stats[m][0]),
                _TIME_FMT.format(stats[m][1]),
            ]
            for m in markers
        ]
        return "\n".join(
            [
                f"## {single_file_code}",
                "",
                _md_table(["Milestone", "Mean (min)", "Std Dev (min)"], table_rows),
            ]
        )
    else:
        headers = ["Run"] + marker_labels
        table_rows = [
            [rc]
            + [
                f"{rows_dict[rc][m][0]:.2f} \u00b1 {rows_dict[rc][m][1]:.2f}"
                for m in markers
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
# Marker helpers
# ---------------------------------------------------------------------------


def _parse_marker(s, param_name="marker"):
    """Parse a percentage string into an integer milestone value.

    :param s: Percentage string, e.g. ``'50'`` or ``'50%'``.
    :type s: str
    :param param_name: Name used in error messages.
    :type param_name: str
    :return: Integer percentage (0–100).
    :rtype: int
    :raises click.BadParameter: If the value is not a valid integer in [0, 100].
    """
    s = s.rstrip("%")
    try:
        val = int(s)
    except ValueError:
        raise click.BadParameter(
            f"Expected an integer percentage (e.g. '50' or '50%'), got: {s!r}",
            param_hint=param_name,
        )
    if not (0 <= val <= 100):
        raise click.BadParameter(
            f"Marker must be between 0 and 100, got: {val}",
            param_hint=param_name,
        )
    return val


def _marker_label(pct, with_units=True):
    """Human-readable label for a degradation milestone.

    :param pct: Integer percentage milestone.
    :type pct: int
    :param with_units: If ``True``, append ``'\\n(min)'`` for Rich multi-line headers.
    :type with_units: bool
    :return: Label string.
    :rtype: str
    """
    label = f"{pct}%"
    return f"{label}\n(min)" if with_units else label


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_run_deg_times(data_root, run_code, markers, console):
    """Load a Run and compute degradation times for each milestone.

    Delegates to
    :func:`~lysis.analysis.degradation.compute_degradation_marker_stats`.

    :param data_root: Directory containing the HDF5 file.
    :type data_root: str
    :param run_code: Run code (HDF5 filename without extension).
    :type run_code: str
    :param markers: List of integer percentage milestones.
    :type markers: list[int]
    :param console: Rich Console for error messages.
    :return: Dict mapping marker percentage → (mean_min, std_min), or ``None``
        on error.
    :rtype: dict[int, tuple] or None
    """
    from lysis.analysis.degradation import compute_degradation_marker_stats
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.macro_params = run.data.macro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None

    try:
        return compute_degradation_marker_stats(run, markers)
    except Exception as e:
        console.print(f"[red]Error computing times for {run_code}:[/red] {e}")
        return None
    finally:
        run.data.close()


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="deg-time")
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
        "Example: --markdown -, --markdown times.md"
    ),
)
@click.option(
    "--add",
    "add_markers",
    multiple=True,
    metavar="PCT",
    help=(
        "Add a degradation milestone to the table (repeatable).  "
        "PCT is an integer percentage (0–100), optionally with a '%' suffix.  "
        "Example: --add 10 --add 90%"
    ),
)
@click.option(
    "--drop",
    "drop_markers",
    multiple=True,
    metavar="PCT",
    help=(
        "Remove a default milestone from the table (repeatable).  "
        "Example: --drop 5 --drop 100%"
    ),
)
@click.pass_context
def deg_time(ctx, path, sort_mode, no_progress, markdown_out, add_markers, drop_markers):
    """Print degradation-time tables for one or more simulation Runs.

    Computes the mean and standard deviation of the time (minutes) to reach
    each degradation milestone across all macroscale simulations.  Mirrors
    the notebook's "## Degradation Time Table".

    PATH may be a single HDF5 file or a directory.  In directory mode every
    .h5 file is processed and results are shown as a table with one row per
    Run and one column per milestone.

    Default milestones: 5%, 20%, 50%, 80%, 100%.

    \b
    Examples:
        lysis deg-time data/TB-xi__1_582_867.h5
        lysis deg-time data/lysis-front/
        lysis deg-time data/ --add 10 --drop 5
        lysis deg-time data/ --markdown -
        lysis deg-time data/ --markdown times.md
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

    # Build effective marker list
    drop_set = set()
    for raw in drop_markers:
        try:
            drop_set.add(_parse_marker(raw))
        except click.BadParameter as e:
            console.print(f"[red]Error:[/red] --drop: {e.format_message()}")
            ctx.exit(1)
            return

    markers = [m for m in _DEFAULT_MARKERS if m not in drop_set]

    for raw in add_markers:
        try:
            m = _parse_marker(raw)
        except click.BadParameter as e:
            console.print(f"[red]Error:[/red] --add: {e.format_message()}")
            ctx.exit(1)
            return
        if m not in markers:
            markers.append(m)

    markers = sorted(markers)

    if not markers:
        console.print("[red]Error:[/red] No milestones remain after applying --drop.")
        ctx.exit(1)
        return

    # -----------------------------------------------------------------------
    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)

        if not no_progress:
            with console.status(f"Computing degradation times for {run_code}..."):
                stats = _load_run_deg_times(data_root, run_code, markers, console)
        else:
            stats = _load_run_deg_times(data_root, run_code, markers, console)

        if stats is None:
            ctx.exit(1)
            return

        if markdown_out is not None:
            _emit_markdown(
                _deg_time_to_markdown({run_code: stats}, markers, run_code),
                markdown_out,
                console,
            )
        else:
            console.print(f"[bold]{run_code}[/bold]\n")
            for m in markers:
                mean, std = stats[m]
                console.print(
                    f"  {m}%: "
                    f"{_TIME_FMT.format(mean)} \u00b1 {_TIME_FMT.format(std)} min"
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
                prog.add_task("Computing degradation times...", total=len(run_codes))
                if prog is not None
                else None
            )
            for run_code in run_codes:
                stats = _load_run_deg_times(path, run_code, markers, console)
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
                _deg_time_to_markdown(rows, markers),
                markdown_out,
                console,
            )
        else:
            table = Table(show_header=True, header_style="bold", show_lines=True)
            table.add_column("Run", style="bold", no_wrap=True)
            for m in markers:
                table.add_column(_marker_label(m), justify="right")

            for run_code in ordered:
                stats = rows[run_code]
                cells = [run_code]
                for m in markers:
                    mean, std = stats[m]
                    cells.append(
                        f"{_TIME_FMT.format(mean)}\n\u00b1 {_TIME_FMT.format(std)}"
                    )
                table.add_row(*cells)

            console.print(table)
