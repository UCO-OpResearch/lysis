"""``lysis compare`` — compare two Runs using the 2-sample KS test.

Thin wrapper around :func:`lysis.analysis.compare.compare_runs_ks`.  Opens
the two HDF5 files named on the command line, dispatches on ``--which``,
and prints the :func:`scipy.stats.ks_2samp` statistic and p-value for each
measure in the selected set.
"""

import os

import click

from lysis.analysis.compare import MEASURE_EXTRACTORS, compare_runs_ks
from lysis.cli import cli


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


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
@click.argument(
    "file1",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
@click.argument(
    "file2",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
@click.pass_context
def compare(ctx, which, file1, file2):
    """Compare two Runs using the 2-sample Kolmogorov-Smirnov test.

    For each measure in the set selected by ``--which``, the per-simulation
    arrays from FILE1 and FILE2 are compared via
    :func:`scipy.stats.ks_2samp` and the statistic and p-value are printed.

    \b
    Examples:
        lysis compare --which micro-stats run_A.h5 run_B.h5
    """
    console = ctx.obj["console"]

    file1 = os.path.abspath(file1)
    file2 = os.path.abspath(file2)

    run_code_1 = os.path.splitext(os.path.basename(file1))[0]
    run_code_2 = os.path.splitext(os.path.basename(file2))[0]

    run1 = _open_run(os.path.dirname(file1), run_code_1, console)
    if run1 is None:
        ctx.exit(1)
        return

    try:
        run2 = _open_run(os.path.dirname(file2), run_code_2, console)
        if run2 is None:
            ctx.exit(1)
            return

        try:
            results = compare_runs_ks(run1, run2, which)
        except Exception as e:
            console.print(f"[red]Error comparing runs:[/red] {e}")
            ctx.exit(1)
            return
        finally:
            run2.data.close()
    finally:
        run1.data.close()

    console.print(
        f"Comparing [bold]{run_code_1}[/bold] vs [bold]{run_code_2}[/bold] "
        f"using [cyan]{which}[/cyan] measures:\n"
    )
    for label, result in results.items():
        console.print(
            f"  {label}: statistic = {result.statistic:.6f}, "
            f"p-value = {result.pvalue:.6g}"
        )
