"""``lysis summarize`` — print summary statistics for a simulation Run."""

import os

import click

from lysis.cli import cli


@cli.command()
@click.argument("hdf5_file", type=click.Path(exists=True, dir_okay=False))
@click.pass_context
def summarize(ctx, hdf5_file):
    """Print summary statistics for a simulation Run.

    Reads the HDF5 file at HDF5_FILE, computes six key statistics across all
    macroscale simulations, and prints each as ``mean ± standard deviation``.

    \b
    Example:
        lysis summarize data/2024-04-16-1900/2024-04-16-1900.h5
    """
    from lysis.analysis.degradation import compute_run_statistics
    from lysis.config.run import Run

    console = ctx.obj["console"]

    hdf5_file = os.path.abspath(hdf5_file)
    os_path = os.path.dirname(hdf5_file)
    run_code = os.path.splitext(os.path.basename(hdf5_file))[0]
    data_root = os.path.dirname(os_path)

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.macro_params = run.data.macro_params
    except Exception as e:
        console.print(f"[red]Error loading run:[/red] {e}")
        ctx.exit(1)
        return

    try:
        stats = compute_run_statistics(run)
    except Exception as e:
        console.print(f"[red]Error computing statistics:[/red] {e}")
        ctx.exit(1)
        return
    finally:
        run.data.close()

    console.print(f"[bold]{run_code}[/bold]\n")
    for metric in stats.index.get_level_values(0).unique():
        x = stats[metric]
        console.print(
            f"  {metric}: {x['Mean']:,.3f} \u00b1 {x['Standard Deviation']:,.3f}"
        )
