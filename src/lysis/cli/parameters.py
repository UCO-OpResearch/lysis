"""``lysis parameters`` — print parameter tables for simulation Runs."""

import os

import click

from lysis.cli import cli

# Each entry: (attr_name, source, display_units, fmt)
# source:        "macro" | "micro" | "computed" — informational only
# display_units: pint unit string to convert to before showing the magnitude,
#                or None to show the value's magnitude in its natural units
# fmt:           Python format string applied to the final numeric value
_DEFAULT_PARAMS = [
    # ---- Macro physical -----------------------------------------------
    ("pore_size",             "macro",    "microns",   "{:.4f}"),
    ("diffusion_coeff",       "macro",    None,         "{:1.5e}"),
    ("forced_unbind",         "macro",    None,         "{:.3f}"),
    ("average_bound_time",    "macro",    None,         "{:.2f}"),
    ("grid_node_distance",    "computed", "microns",    "{:.3f}"),
    # ---- Grid geometry ------------------------------------------------
    ("cols",                  "macro",    None,         "{:,}"),
    ("rows",                  "macro",    None,         "{:,}"),
    ("fiber_rows",            "macro",    None,         "{:,}"),
    ("empty_rows",            "macro",    None,         "{:,}"),
    ("full_row",              "macro",    None,         "{:,}"),
    ("xz_row",                "macro",    None,         "{:,}"),
    ("total_edges",           "macro",    None,         "{:,}"),
    ("total_fibers",          "macro",    None,         "{:,}"),
    ("total_molecules",       "macro",    None,         "{:,}"),
    ("moving_probability",    "macro",    None,         "{:.2f}"),
    # ---- Macro run controls -------------------------------------------
    ("macro_simulations",     "macro",    None,         "{:,}"),
    ("total_time",            "macro",    "minutes",    "{:.0f}"),
    ("time_step",             "macro",    None,         "{:1.3e}"),
    ("total_time_steps",      "macro",    None,         "{:,}"),
    ("macro_seed",            "macro",    None,         "{:,}"),
    ("save_interval",         "macro",    None,         "{:.0f}"),
    ("number_of_saves",       "macro",    None,         "{:,}"),
    # ---- Micro physical -----------------------------------------------
    ("bind_rate_tPA",         "micro",    None,         "{:.3f}"),
    ("bind_rate_PLG",         "micro",    None,         "{:.3f}"),
    ("conc_free_PLG",         "micro",    None,         "{:.2f}"),
    ("diss_const_tPA_wPLG",  "micro",    None,         "{:.4f}"),
    ("diss_const_tPA_woPLG", "micro",    None,         "{:.4f}"),
    ("diss_const_PLG_intact", "micro",   None,         "{:.2f}"),
    ("diss_const_PLG_nicked", "micro",   None,         "{:.2f}"),
    ("deg_rate_fibrin",       "micro",   None,          "{:.2f}"),
    ("fiber_radius",          "micro",   "nanometers",  "{:.4f}"),
    ("binding_sites",         "micro",   None,          "{:.6f}"),
    # ---- Micro run controls -------------------------------------------
    ("micro_simulations",     "micro",   None,          "{:,}"),
    ("micro_seed",            "micro",   None,          "{:,}"),
]

# Computed (virtual) parameters not directly stored on any params object
_COMPUTED = {
    "grid_node_distance": lambda mp: (
        mp.pore_size + 2 * mp.micro_params.fiber_radius
    ),
}


def _get_raw(macro_params, attr_name):
    """Return the raw value of a parameter (Pint Quantity or scalar), or None."""
    if attr_name in _COMPUTED:
        return _COMPUTED[attr_name](macro_params)
    if hasattr(macro_params, attr_name):
        return getattr(macro_params, attr_name)
    micro = getattr(macro_params, "micro_params", None)
    if micro is not None and hasattr(micro, attr_name):
        return getattr(micro, attr_name)
    return None


def _format_value(raw, display_units, fmt):
    """Convert and format a raw parameter value for table display."""
    if raw is None:
        return "N/A"
    if hasattr(raw, "to"):  # Pint Quantity
        num = raw.to(display_units).magnitude if display_units else raw.magnitude
        return fmt.format(num) if fmt else str(num)
    return fmt.format(raw) if fmt else str(raw)


def _auto_format(raw):
    """Format a value for --add parameters (no explicit format or units hint)."""
    if raw is None:
        return "N/A"
    if hasattr(raw, "to"):  # Pint Quantity — show with natural units
        try:
            return f"{raw:.4g~P}"
        except Exception:
            return str(raw)
    if isinstance(raw, int):
        return f"{raw:,}"
    if isinstance(raw, float):
        return f"{raw:.4g}"
    return str(raw)


def _param_header(attr_name, display_units, natural_units):
    """Build a single-line row label with units in parentheses."""
    units = display_units or natural_units.get(attr_name)
    return f"{attr_name} ({units})" if units else attr_name


def _load_run_params(data_root, run_code, param_specs, add_names, console):
    """Load a Run and return a dict of formatted parameter values, or None on error."""
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        run.macro_params = run.data.macro_params
        macro_params = run.macro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None

    try:
        values = {}
        for attr_name, _src, display_units, fmt in param_specs:
            raw = _get_raw(macro_params, attr_name)
            values[attr_name] = _format_value(raw, display_units, fmt)
        for attr_name in add_names:
            raw = _get_raw(macro_params, attr_name)
            values[attr_name] = _auto_format(raw)
        return values
    except Exception as e:
        console.print(f"[red]Error reading parameters for {run_code}:[/red] {e}")
        return None
    finally:
        run.data.close()


@cli.command()
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
    "--add",
    "add_params",
    multiple=True,
    metavar="PARAM",
    help=(
        "Add a parameter to the table (repeatable). "
        "PARAM must be a MacroParameters or MicroParameters attribute name. "
        "Example: --add protofibril_radius --add activation_rate_PLG"
    ),
)
@click.option(
    "--drop",
    "drop_params",
    multiple=True,
    metavar="PARAM",
    help=(
        "Remove a parameter from the default table (repeatable). "
        "Example: --drop empty_rows --drop macro_seed"
    ),
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
        "Example: --markdown -, --markdown params.md"
    ),
)
@click.pass_context
def parameters(ctx, path, sort_mode, no_progress, add_params, drop_params, markdown_out):
    """Print parameter tables for one or more simulation Runs.

    PATH may be a single HDF5 file or a directory. When a directory is given,
    parameters for every .h5 file are printed as a table with one column per
    Run and one row per parameter. Parameters are grouped into macro and micro
    sections, separated by a divider.

    \b
    Examples:
        lysis parameters data/TB-xi__1_582_867.h5
        lysis parameters data/lysis-front/
        lysis parameters data/lysis-front/ --sort alpha
        lysis parameters data/ --drop empty_rows --add protofibril_radius
        lysis parameters data/ --markdown -
        lysis parameters data/ --markdown params.md
    """
    from contextlib import nullcontext

    from lysis.analysis.summary import parameters_table
    from lysis.tools.display import emit_markdown, params_df_to_markdown, params_df_to_rich
    from lysis.config.parameters import MacroParameters

    console = ctx.obj["console"]
    path = os.path.abspath(path)

    # Suppress progress when producing structured markdown output
    if markdown_out is not None:
        no_progress = True

    # Build effective parameter spec list
    drop_set = set(drop_params)
    param_specs = [p for p in _DEFAULT_PARAMS if p[0] not in drop_set]
    add_names = list(add_params)

    # Pre-build natural-units dict for column headers (no run needed)
    natural_units = MacroParameters.units()

    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)

        if not no_progress:
            with console.status(f"Loading parameters for {run_code}..."):
                values = _load_run_params(
                    data_root, run_code, param_specs, add_names, console
                )
        else:
            values = _load_run_params(
                data_root, run_code, param_specs, add_names, console
            )

        if values is None:
            ctx.exit(1)
            return

        df = parameters_table(
            {run_code: values}, param_specs, add_names, natural_units, [run_code]
        )

        if markdown_out is not None:
            emit_markdown(
                params_df_to_markdown(df, single_code=run_code),
                markdown_out,
                console,
            )
        else:
            console.print(params_df_to_rich(df))

    else:
        # --- directory mode: parameters × runs table ---
        from lysis.tools.runcode_sort import smart_sort
        from rich.progress import (
            BarColumn,
            MofNCompleteColumn,
            Progress,
            TextColumn,
            TimeRemainingColumn,
        )

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
                prog.add_task("Loading parameters...", total=len(run_codes))
                if prog is not None
                else None
            )
            for run_code in run_codes:
                vals = _load_run_params(
                    path, run_code, param_specs, add_names, console
                )
                if prog is not None:
                    prog.advance(_task)
                if vals is not None:
                    rows[run_code] = vals

        if not rows:
            ctx.exit(1)
            return

        # Preserve sort order, skipping any failed runs
        ordered = [rc for rc in run_codes if rc in rows]
        df = parameters_table(rows, param_specs, add_names, natural_units, ordered)

        if markdown_out is not None:
            emit_markdown(
                params_df_to_markdown(df),
                markdown_out,
                console,
            )
        else:
            console.print(params_df_to_rich(df))
