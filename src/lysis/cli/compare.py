"""``lysis compare`` — compare Runs across two folders using the 2-sample KS test.

Thin wrapper around :func:`lysis.analysis.compare.compare_runs`.  Takes
two folders, finds the set of run codes present in both, and reports the
:func:`scipy.stats.ks_2samp` statistic and p-value for each measure in the
set selected by ``--which`` plus a symmetric percent-difference on each
paired scalar summary.
"""

import os

import click

from lysis.analysis.compare import MEASURE_EXTRACTORS, compare_runs
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
    :param which: Key into :data:`MEASURE_EXTRACTORS`.
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
    "--which",
    type=click.Choice(list(MEASURE_EXTRACTORS.keys()), case_sensitive=False),
    required=True,
    help="Which set of measures to compare.",
)
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
@click.argument(
    "folder1",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.argument(
    "folder2",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.pass_context
def compare(ctx, which, sort_mode, no_progress, markdown_out, folder1, folder2):
    """Compare Runs across two folders using the 2-sample Kolmogorov-Smirnov test.

    For every run code whose ``.h5`` file appears in both FOLDER1 and
    FOLDER2 the command computes, per measure set selected by ``--which``:

    - a 2-sample :func:`scipy.stats.ks_2samp` test on each per-simulation
      array pair (``KS`` and ``p-value`` columns);
    - a symmetric percent difference on each scalar summary stat (``% diff``
      columns, signed so positive means FOLDER2 > FOLDER1).

    Results are rendered as a table with one row per Run.

    \b
    Examples:
        lysis compare --which micro-stats data/runA/ data/runB/
        lysis compare --which macro-stats data/runA/ data/runB/
        lysis compare --which macro-stats data/runA/ data/runB/ --sort alpha
        lysis compare --which macro-stats data/runA/ data/runB/ --no-progress
        lysis compare --which macro-stats data/runA/ data/runB/ --markdown -
        lysis compare --which macro-stats data/runA/ data/runB/ --markdown out.md
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import compare_stats_table
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
