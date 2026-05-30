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

from .._metadata import __version__


# Shared console — auto-detects TTY vs file redirection.
# In a terminal: colors and formatting.  Redirected to file: plain text.
console = Console()


@click.group()
@click.version_option(__version__, prog_name="lysis")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv).")
@click.pass_context
def cli(ctx, verbose):
    """Lysis simulation data tools."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["console"] = console

    # Resolve src/lysis/ version + dirty state once per invocation.
    # Cached on ctx.obj so individual subcommands consult the cache
    # instead of reshelling to git.
    from lysis.tools.provenance.execution import (  # noqa: PLC0415
        _resolve_dirty,
        _resolve_version,
    )
    from lysis.tools.provenance import mark_dirty_warning_emitted  # noqa: PLC0415

    ctx.obj["lysis_version"] = _resolve_version()
    ctx.obj["lysis_dirty"] = _resolve_dirty()  # "clean"|"dirty"|"unknown"

    if ctx.obj["lysis_dirty"] == "dirty":
        # Route to stderr so it doesn't pollute the stdout output of
        # any subcommand (e.g. machine-readable markdown tables).  Use
        # a dedicated stderr Console rather than ``console`` (which is
        # stdout-attached) so colours are correctly disabled when
        # stderr is a pipe / file.
        from rich.console import Console as _Console  # noqa: PLC0415
        _Console(stderr=True).print(
            "[yellow]Warning:[/yellow] src/lysis/ has uncommitted changes "
            f"(version stamped as {ctx.obj['lysis_version']}).",
            highlight=False,
        )
        # Suppress the redundant inline warnings.warn() inside subsequent
        # gather_pipeline_provenance() / gather_init_provenance() calls
        # in this process.  Non-CLI callers (notebooks, scripts) never
        # call this, so they still receive the warning naturally.
        mark_dirty_warning_emitted()


def main():
    """Entry point for the ``lysis`` console command."""
    cli()


# Import command modules to register them with the Click group.
from lysis.cli import compare, convert, deg_rate, deg_time, diff, init_experiment, init_macroscale, macro_stats, micro_stats, parameters, rename, run_macro, run_micro, validate  # noqa: E402, F401
