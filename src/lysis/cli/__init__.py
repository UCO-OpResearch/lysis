"""Lysis command-line interface.

Usage::

    lysis <command> [options]

Available commands::

    convert     Convert simulation data between specification formats
    validate    Validate simulation data against a specification

Use ``lysis --help`` or ``lysis <command> --help`` for details.
"""

import click
from rich.console import Console


# Shared console — auto-detects TTY vs file redirection.
# In a terminal: colors and formatting.  Redirected to file: plain text.
console = Console()


@click.group()
@click.version_option("0.1.0", prog_name="lysis")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv).")
@click.pass_context
def cli(ctx, verbose):
    """Lysis simulation data tools."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["console"] = console


def main():
    """Entry point for the ``lysis`` console command."""
    cli()


# Import command modules to register them with the Click group.
from lysis.cli import compare, convert, deg_rate, deg_time, diff, init_experiment, init_macroscale, macro_stats, micro_stats, parameters, run_macro, run_micro, validate  # noqa: E402, F401
