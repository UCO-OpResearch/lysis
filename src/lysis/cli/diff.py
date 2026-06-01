"""``lysis diff`` — element-wise comparison of HDF5 data tables.

Thin wrapper around :func:`lysis.analysis.diff.diff_runs`.  Takes two
paths (folders or ``.h5`` files) and compares the on-disk microscale
and/or macroscale datasets element-wise with a 2-ULP floating-point
tolerance, reporting per-table match / diff / shape-mismatch /
missing.  Which scales are compared is auto-detected from the two
Runs' data collections: if either Run is missing a scale the other
has, that scale is skipped and a yellow warning is printed; if the two
Runs share no scales at all the command errors out.

In file-pair mode ``--table TABLE`` switches to a side-by-side
page-through view of a single dataset, with differing rows
highlighted.  This option has no effect in folder-pair mode.
"""

import os

import click

from lysis.analysis.diff import diff_runs, list_tables
from lysis.cli import cli
from lysis.cli.compare import (
    _list_h5_run_codes,
    _open_run,
    _split_h5_path,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _diff_one(folder1, folder2, run_code, console):
    """Open both Runs for *run_code*, call :func:`diff_runs`, then close.

    :param folder1: Directory containing the first Run's HDF5 file.
    :type folder1: str
    :param folder2: Directory containing the second Run's HDF5 file.
    :type folder2: str
    :param run_code: Run code shared by both folders.
    :type run_code: str
    :param console: Rich Console for error / warning messages.
    :return: Result dict from :func:`diff_runs`, or ``None`` on error.
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
            return diff_runs(run1, run2)
        except ValueError as e:
            console.print(f"[red]Error diffing {run_code}:[/red] {e}")
            return None
        except Exception as e:
            console.print(f"[red]Error diffing {run_code}:[/red] {e}")
            return None
        finally:
            run2.data.close()
    finally:
        run1.data.close()


def _diff_file_pair(path1, path2, console):
    """Open two ``.h5`` files and call :func:`diff_runs` on them.

    :param path1: Path to the first HDF5 file.
    :type path1: str
    :param path2: Path to the second HDF5 file.
    :type path2: str
    :param console: Rich Console for error messages.
    :return: Result dict from :func:`diff_runs`, or ``None`` on error.
    :rtype: dict or None
    """
    root1, code1 = _split_h5_path(path1)
    root2, code2 = _split_h5_path(path2)
    run1 = _open_run(root1, code1, console)
    if run1 is None:
        return None
    try:
        run2 = _open_run(root2, code2, console)
        if run2 is None:
            return None
        try:
            return diff_runs(run1, run2)
        except ValueError as e:
            console.print(f"[red]Error:[/red] {e}")
            return None
        except Exception as e:
            console.print(f"[red]Error diffing files:[/red] {e}")
            return None
        finally:
            run2.data.close()
    finally:
        run1.data.close()


def _render_side_by_side_from_files(path1, path2, table_label, console, ctx):
    """Open two ``.h5`` files and render the named dataset side-by-side.

    Uses the scales present in both Runs (via
    :func:`~lysis.analysis.diff.available_scales`) and looks up
    *table_label* in each Run.  Errors out if the label is not present
    in either Run or if the two runs cannot be opened.

    :param path1: Path to the first HDF5 file.
    :type path1: str
    :param path2: Path to the second HDF5 file.
    :type path2: str
    :param table_label: Dataset label to render.
    :type table_label: str
    :param console: Rich Console.
    :param ctx: Click context (used for ``ctx.exit`` on failure).
    """
    from lysis.tools.display import render_dataset_side_by_side

    root1, code1 = _split_h5_path(path1)
    root2, code2 = _split_h5_path(path2)
    run1 = _open_run(root1, code1, console)
    if run1 is None:
        ctx.exit(1)
        return
    try:
        run2 = _open_run(root2, code2, console)
        if run2 is None:
            ctx.exit(1)
            return
        try:
            labels1 = list_tables(run1)
            labels2 = list_tables(run2)

            if table_label not in labels1 or table_label not in labels2:
                available = sorted(set(labels1) | set(labels2))
                missing_from = []
                if table_label not in labels1:
                    missing_from.append(os.path.basename(path1))
                if table_label not in labels2:
                    missing_from.append(os.path.basename(path2))
                console.print(
                    f"[red]Error:[/red] table {table_label!r} not found in "
                    f"{', '.join(missing_from)}."
                )
                console.print(
                    "Available tables:\n  " + "\n  ".join(available)
                )
                ctx.exit(1)
                return

            arr1 = _read_table(run1, table_label)
            arr2 = _read_table(run2, table_label)
            render_dataset_side_by_side(
                arr1,
                arr2,
                table_label,
                console,
                file1_name=os.path.basename(path1),
                file2_name=os.path.basename(path2),
            )
        finally:
            run2.data.close()
    finally:
        run1.data.close()


def _read_table(run, label):
    """Read a single dataset out of *run* by display label.

    Accepts labels of the form ``microscale_out/<name>`` or
    ``macroscale_out[<sim>]/<name>`` — the same format produced by
    :func:`~lysis.analysis.diff.list_tables`.

    :param run: Run with data open.
    :type run: Run
    :param label: Display label.
    :type label: str
    :return: Array copy of the dataset.
    :rtype: numpy.ndarray
    :raises KeyError: If *label* is not in the Run's table list.
    """
    data = run.data
    if label.startswith("microscale_out/"):
        name = label[len("microscale_out/"):]
        arr = getattr(data.microscale_out, name)
    elif label.startswith("macroscale_out["):
        idx_end = label.index("]")
        sim = int(label[len("macroscale_out["):idx_end])
        name = label[idx_end + 2:]  # skip "]/"
        arr = getattr(data.macroscale_out[sim], name)
    else:
        raise KeyError(label)
    if arr is None:  # absent optional dataset -- not a comparable table
        raise KeyError(label)
    return arr[:]


# ---------------------------------------------------------------------------
# Click command
# ---------------------------------------------------------------------------


@cli.command(name="diff")
@click.option(
    "--sort",
    "sort_mode",
    type=click.Choice(["smart", "alpha"], case_sensitive=False),
    default="smart",
    show_default=True,
    help=(
        "Sort order for the common run codes in folder-pair mode. "
        "'smart' detects roman numerals and python-style integers and "
        "sorts them numerically; 'alpha' uses plain lexicographic order."
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
        "or provide a filename to write to a file.  Progress indicators "
        "are suppressed automatically.  Example: --markdown -, "
        "--markdown diff.md"
    ),
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help=(
        "Folder-pair mode: emit the full per-table result matrix "
        "instead of the compact summary (one row per Run).  No effect "
        "in file-pair mode."
    ),
)
@click.option(
    "--table",
    "table_label",
    type=str,
    default=None,
    metavar="TABLE",
    help=(
        "File-pair mode only: render TABLE from both files side-by-side "
        "in a pager, with differing rows highlighted.  TABLE is a "
        "dataset label such as 'microscale_out/tpa_leaving_time' or "
        "'macroscale_out[00]/snapshot_time'.  Skips the full diff; the "
        "named dataset is rendered instead."
    ),
)
@click.argument(
    "path1",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)
@click.argument(
    "path2",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)
@click.pass_context
def diff(
    ctx,
    sort_mode,
    no_progress,
    markdown_out,
    verbose,
    table_label,
    path1,
    path2,
):
    """Diff HDF5 data tables between two Runs or two folders of Runs.

    PATH1 and PATH2 may both be directories (folder-pair mode) or both
    be ``.h5`` files (file-pair mode); mixed input is rejected.

    Which scales (``microscale_out`` / ``macroscale_out``) get compared
    is auto-detected from each Run's data collections.  If one Run
    lacks a scale the other has, that scale is skipped and a yellow
    warning is printed; if neither Run has any scale in common, the
    command errors out.

    \b
    Folder-pair mode (the default use case):
        For every run code whose ``.h5`` file appears in both PATH1 and
        PATH2 the command diffs the two Runs.  By default emits a
        compact summary per Run (``Result`` / ``Max % Diff`` /
        ``Worst Table``); ``--verbose`` switches to the full per-table
        matrix.

    \b
    File-pair mode:
        PATH1 and PATH2 are each diffed as a single Run, producing a
        detailed per-table breakdown (one row per dataset with
        ``Status``, ``Mismatches``, ``Max % Diff``, ``Location``).
        ``--table TABLE`` switches to a side-by-side page-through view
        of a single dataset (differing rows highlighted).  ``--table``
        is file-pair only.

    \b
    Examples:
        lysis diff data/runA/ data/runB/
        lysis diff data/runA/ data/runB/ --verbose --markdown out.md
        lysis diff runA.h5 runB.h5
        lysis diff runA.h5 runB.h5 --markdown out.md
        lysis diff runA.h5 runB.h5 --table microscale_out/tpa_leaving_time
        lysis diff runA.h5 runB.h5 --table "macroscale_out[00]/snapshot_time"
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import diff_detail_table, diff_summary_table
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
    path1 = os.path.abspath(path1)
    path2 = os.path.abspath(path2)

    if markdown_out is not None:
        no_progress = True

    path1_is_file = os.path.isfile(path1)
    path2_is_file = os.path.isfile(path2)
    if path1_is_file != path2_is_file:
        console.print(
            "[red]Error:[/red] PATH1 and PATH2 must both be directories or "
            "both be .h5 files."
        )
        ctx.exit(1)
        return

    if table_label is not None and not (path1_is_file and path2_is_file):
        console.print(
            "[red]Error:[/red] --table is only supported in file-pair mode "
            "(pass two .h5 files, not directories)."
        )
        ctx.exit(1)
        return

    if path1_is_file and path2_is_file:
        if not (path1.lower().endswith(".h5") and path2.lower().endswith(".h5")):
            console.print(
                "[red]Error:[/red] file-pair inputs must have .h5 extensions."
            )
            ctx.exit(1)
            return

        if table_label is not None:
            _render_side_by_side_from_files(
                path1, path2, table_label, console, ctx
            )
            return

        result = _diff_file_pair(path1, path2, console)
        if result is None:
            ctx.exit(1)
            return

        for missing in result.get("scales_skipped_in_run1", []):
            console.print(
                f"[yellow]Warning:[/yellow] {missing} is not present in "
                f"{os.path.basename(path1)}; comparing only the shared scales."
            )
        for missing in result.get("scales_skipped_in_run2", []):
            console.print(
                f"[yellow]Warning:[/yellow] {missing} is not present in "
                f"{os.path.basename(path2)}; comparing only the shared scales."
            )

        df = diff_detail_table(result["data_diff"])
        if markdown_out is not None:
            emit_markdown(stats_df_to_markdown(df, "Table"), markdown_out, console)
        else:
            console.print(stats_df_to_rich(df, "Table"))
        return

    codes1 = set(_list_h5_run_codes(path1))
    codes2 = set(_list_h5_run_codes(path2))
    common = codes1 & codes2

    if not common:
        console.print(
            f"[yellow]No run codes in common between {path1} and {path2}[/yellow]"
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
    skipped_in_path1 = set()
    skipped_in_path2 = set()
    with _pctx as prog:
        _task = (
            prog.add_task("Diffing runs...", total=len(common))
            if prog is not None
            else None
        )
        for run_code in common:
            try:
                result = _diff_one(path1, path2, run_code, console)
                if result is not None:
                    all_results[run_code] = result
                    skipped_in_path1.update(result.get("scales_skipped_in_run1", []))
                    skipped_in_path2.update(result.get("scales_skipped_in_run2", []))
            finally:
                if prog is not None:
                    prog.advance(_task)

    if not all_results:
        ctx.exit(1)
        return

    ordered = {rc: all_results[rc] for rc in common if rc in all_results}
    if verbose:
        from lysis.analysis.summary import compare_stats_table

        df = compare_stats_table(ordered)
    else:
        df = diff_summary_table(ordered)

    if markdown_out is not None:
        emit_markdown(stats_df_to_markdown(df, "Run"), markdown_out, console)
    else:
        console.print(stats_df_to_rich(df, "Run"))

    for missing in sorted(skipped_in_path1):
        console.print(
            f"[yellow]Warning:[/yellow] {missing} was missing from one or "
            f"more Runs in {os.path.basename(path1)}; those scales were "
            "skipped."
        )
    for missing in sorted(skipped_in_path2):
        console.print(
            f"[yellow]Warning:[/yellow] {missing} was missing from one or "
            f"more Runs in {os.path.basename(path2)}; those scales were "
            "skipped."
        )
