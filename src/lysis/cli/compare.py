"""``lysis compare`` — compare Runs across two folders using the 2-sample KS test.

Thin wrapper around :func:`lysis.analysis.compare.compare_runs`.  Takes
a measure-set name and two folders, finds the set of run codes present
in both, and reports the :func:`scipy.stats.ks_2samp` statistic and
p-value for each measure in the selected set plus a symmetric
percent-difference on each paired scalar summary.
"""

import os

import click

from lysis.analysis.compare import (
    DATA_TABLE_EXTRACTORS,
    available_measure_sets,
    compare_runs,
)
from lysis.cli import cli


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _list_h5_run_codes(folder):
    """Return the run codes for every ``.h5`` file directly inside *folder*.

    :param folder: Directory to scan.
    :type folder: str
    :return: Run codes (HDF5 filenames without extension).
    :rtype: list[str]
    """
    return [
        os.path.splitext(f)[0]
        for f in os.listdir(folder)
        if f.lower().endswith(".h5")
    ]


def _open_run(data_root, run_code, console):
    """Open a Run for reading, or return ``None`` on error.

    :param data_root: Directory containing the HDF5 file.
    :type data_root: str
    :param run_code: Run code (HDF5 filename without extension).
    :type run_code: str
    :param console: Rich Console for error messages.
    :return: Opened Run, or ``None`` on failure.
    :rtype: lysis.config.run.Run or None
    """
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.micro_params = run.data.micro_params
        if run.data.macro_params is not None:
            run.macro_params = run.data.macro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None
    return run


def _compare_one(folder1, folder2, run_code, which, console):
    """Open both Runs for *run_code*, run :func:`compare_runs`, then close.

    :param folder1: Directory containing the first Run's HDF5 file.
    :type folder1: str
    :param folder2: Directory containing the second Run's HDF5 file.
    :type folder2: str
    :param run_code: Run code shared by both folders.
    :type run_code: str
    :param which: Measure-set key (see :func:`available_measure_sets`).
    :type which: str
    :param console: Rich Console for error messages.
    :return: Result dict from :func:`compare_runs`, or ``None`` on error.
    :rtype: dict or None
    """
    run1 = _open_run(folder1, run_code, console)
    if run1 is None:
        return None
    try:
        run2 = _open_run(folder2, run_code, console)
        if run2 is None:
            return None
        try:
            return compare_runs(run1, run2, which)
        except Exception as e:
            console.print(f"[red]Error comparing {run_code}:[/red] {e}")
            return None
        finally:
            run2.data.close()
    finally:
        run1.data.close()


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="compare")
@click.option(
    "--sort",
    "sort_mode",
    type=click.Choice(["smart", "alpha"], case_sensitive=False),
    default="smart",
    show_default=True,
    help=(
        "Sort order for the common run codes. "
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
        "Example: --markdown -, --markdown compare.md"
    ),
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help=(
        "For 'micro-data' and 'data' comparisons, emit the full per-table "
        "result matrix instead of the compact summary (one row per Run).  "
        "No effect on 'micro-stats' / 'macro-stats' modes."
    ),
)
@click.argument(
    "which",
    type=click.Choice(available_measure_sets(), case_sensitive=False),
)
@click.argument(
    "folder1",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.argument(
    "folder2",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.pass_context
def compare(
    ctx, sort_mode, no_progress, markdown_out, verbose, which, folder1, folder2
):
    """Compare Runs across two folders.

    For every run code whose ``.h5`` file appears in both FOLDER1 and
    FOLDER2 the command compares the two Runs.  The comparison performed
    depends on WHICH:

    \b
    - ``micro-stats`` / ``macro-stats`` — 2-sample Kolmogorov-Smirnov
      tests on per-simulation arrays plus symmetric percent differences
      on scalar summary stats (signed so positive means FOLDER2 > FOLDER1).
    - ``micro-data`` / ``data`` — element-wise exact-match check on every
      non-log data table in the paired HDF5 files.  HDF5 attributes and
      log tables (``micro_log`` / ``macro_log``) are never examined.
      By default emits a compact summary per Run: ``Exact Match`` if
      every table matched, or the count of differing tables plus the
      maximum symmetric percent difference and the dataset label of
      the worst mismatch.  Use ``--verbose`` to emit a column per
      dataset instead.

    Results are rendered as a table with one row per Run.

    \b
    Examples:
        lysis compare micro-stats data/runA/ data/runB/
        lysis compare macro-stats data/runA/ data/runB/
        lysis compare micro-data data/runA/ data/runB/
        lysis compare data data/runA/ data/runB/
        lysis compare data data/runA/ data/runB/ --verbose --markdown out.md
        lysis compare macro-stats data/runA/ data/runB/ --sort alpha
        lysis compare macro-stats data/runA/ data/runB/ --no-progress
        lysis compare macro-stats data/runA/ data/runB/ --markdown -
        lysis compare macro-stats data/runA/ data/runB/ --markdown out.md
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import (
        compare_data_diff_summary_table,
        compare_stats_table,
    )
    from lysis.tools.display import (
        emit_markdown,
        stats_df_to_markdown,
        stats_df_to_rich,
    )
    from lysis.tools.runcode_sort import smart_sort
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
    )

    console = ctx.obj["console"]
    folder1 = os.path.abspath(folder1)
    folder2 = os.path.abspath(folder2)

    # Suppress progress when producing structured markdown output
    if markdown_out is not None:
        no_progress = True

    codes1 = set(_list_h5_run_codes(folder1))
    codes2 = set(_list_h5_run_codes(folder2))
    common = codes1 & codes2

    if not common:
        console.print(
            f"[yellow]No run codes in common between {folder1} and {folder2}[/yellow]"
        )
        ctx.exit(1)
        return

    common = list(common)
    if sort_mode == "smart":
        common = smart_sort(common)
    else:
        common = sorted(common)

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

    all_results = {}
    with _pctx as prog:
        _task = (
            prog.add_task("Comparing runs...", total=len(common))
            if prog is not None
            else None
        )
        for run_code in common:
            try:
                result = _compare_one(folder1, folder2, run_code, which, console)
                if result is not None:
                    all_results[run_code] = result
            finally:
                if prog is not None:
                    prog.advance(_task)

    if not all_results:
        ctx.exit(1)
        return

    # Preserve sort order, skipping any failed runs
    ordered = {rc: all_results[rc] for rc in common if rc in all_results}
    is_data_diff = which in DATA_TABLE_EXTRACTORS
    if is_data_diff and not verbose:
        df = compare_data_diff_summary_table(ordered)
    else:
        df = compare_stats_table(ordered)

    # Build Rich short-header abbreviations: split the suffix onto a second
    # line so wide measure labels don't make the table unreadable.
    short_headers = {}
    for col in df.columns:
        for suffix in (" KS", " p-value", " % diff"):
            if col.endswith(suffix):
                short_headers[col] = col[: -len(suffix)] + "\n" + suffix.lstrip()
                break

    if markdown_out is not None:
        emit_markdown(stats_df_to_markdown(df, "Run"), markdown_out, console)
    else:
        console.print(stats_df_to_rich(df, "Run", short_headers))
