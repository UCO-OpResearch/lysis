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
    ("unbind_rate_PLi",       "micro",   None,          "{:.2f}"),
    ("activation_rate_PLG",   "micro",   None,          "{:.3f}"),
    ("exposure_rate_binding_site", "micro", None,       "{:.2f}"),
    ("fiber_radius",          "micro",   "nanometers",  "{:.4f}"),
    ("binding_sites",         "micro",   None,          "{:.6f}"),
    ("nodes_in_micro_row",    "micro",   None,          "{:,}"),
    ("snap_proportion",       "micro",   None,          "{:.4f}"),
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


def _get_raw(macro_params, micro_params, attr_name):
    """Return the raw value of a parameter (Pint Quantity or scalar), or None.

    Either ``macro_params`` or ``micro_params`` may be ``None`` when the
    corresponding data collection is absent from the file.  Macroscale and
    computed parameters resolve to ``None`` when no macroscale data is present.
    """
    if attr_name in _COMPUTED:
        return _COMPUTED[attr_name](macro_params) if macro_params is not None else None
    if macro_params is not None and hasattr(macro_params, attr_name):
        return getattr(macro_params, attr_name)
    if micro_params is not None and hasattr(micro_params, attr_name):
        return getattr(micro_params, attr_name)
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


def _default_formatted(param_specs, add_names):
    """Formatted values for freshly-constructed default parameter objects.

    Renders default :class:`~lysis.config.parameters.MicroParameters` and
    :class:`~lysis.config.parameters.MacroParameters` through the same
    formatting path as :func:`_load_run_params`, so a run's formatted value can
    be compared against the model default to decide whether to highlight it.

    :param param_specs: Effective ``(attr_name, source, units, fmt)`` specs.
    :param add_names: Extra attribute names added via ``--add``.
    :return: ``{attr_name: formatted_str}`` for every spec row and add name.
    :rtype: dict[str, str]
    """
    import warnings

    from lysis.config.parameters import MacroParameters, MicroParameters

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        micro = MicroParameters()
        macro = MacroParameters(micro_params=micro)

    defaults = {}
    for attr_name, _src, display_units, fmt in param_specs:
        defaults[attr_name] = _format_value(
            _get_raw(macro, micro, attr_name), display_units, fmt
        )
    for attr_name in add_names:
        defaults[attr_name] = _auto_format(_get_raw(macro, micro, attr_name))
    return defaults


def _nondefault_flags(values, defaults):
    """Map ``{attr_name: bool}`` marking values that differ from the default.

    A value is non-default when its formatted string differs from the default's
    and is not ``"N/A"`` (absent macroscale data renders as ``"N/A"`` and is
    never highlighted).

    :param values: ``{attr_name: formatted_str}`` for one run.
    :param defaults: ``{attr_name: formatted_str}`` from :func:`_default_formatted`.
    :rtype: dict[str, bool]
    """
    return {
        attr: val != "N/A" and val != defaults.get(attr)
        for attr, val in values.items()
    }


def _param_header(attr_name, display_units, natural_units):
    """Build a single-line row label with units in parentheses."""
    units = display_units or natural_units.get(attr_name)
    return f"{attr_name} ({units})" if units else attr_name


def _load_run_params(data_root, run_code, param_specs, add_names, console):
    """Load a Run and read its formatted parameter values.

    :return: ``(values, has_macro)`` where ``values`` is a
        ``{attr_name: formatted_str}`` dict (or ``None`` on error) and
        ``has_macro`` is ``True`` when the file contains a macroscale data
        collection.
    :rtype: tuple[dict | None, bool]
    """
    from lysis.config.run import Run

    try:
        run = Run(data_root, run_code)
        run.open_data()
        # Detect which data collections are present before reading parameters:
        # a microscale-only file has no macroscale parameters (issue #22).
        macro_params = run.data.macro_params
        micro_params = run.data.micro_params
    except Exception as e:
        console.print(f"[red]Error loading {run_code}:[/red] {e}")
        return None, False

    has_macro = macro_params is not None

    if macro_params is None and micro_params is None:
        console.print(f"[yellow]No parameters found for {run_code}[/yellow]")
        run.data.close()
        return None, False

    try:
        values = {}
        for attr_name, _src, display_units, fmt in param_specs:
            raw = _get_raw(macro_params, micro_params, attr_name)
            values[attr_name] = _format_value(raw, display_units, fmt)
        for attr_name in add_names:
            raw = _get_raw(macro_params, micro_params, attr_name)
            values[attr_name] = _auto_format(raw)
        return values, has_macro
    except Exception as e:
        console.print(f"[red]Error reading parameters for {run_code}:[/red] {e}")
        return None, has_macro
    finally:
        run.data.close()


def _macro_param_names():
    """Canonical set of macroscale-only parameter attribute names.

    Derived from the dataclass definitions: a parameter defined on
    :class:`~lysis.config.parameters.MacroParameters` but not on
    :class:`~lysis.config.parameters.MicroParameters` is macroscale-only.  The
    nested ``micro_params`` container is excluded.

    :rtype: set[str]
    """
    import dataclasses

    from lysis.config.parameters import MacroParameters, MicroParameters

    macro = {f.name for f in dataclasses.fields(MacroParameters)}
    micro = {f.name for f in dataclasses.fields(MicroParameters)}
    return (macro - micro) - {"micro_params"}


def _drop_macro_specs(param_specs):
    """Remove macroscale and computed parameter rows from the spec list.

    Used when no run in the table has a macroscale data collection (issue #22),
    so those rows could never carry data.  The macroscale parameters are
    identified from the canonical :class:`MacroParameters` dataclass rather than
    inferred from missing values, so a *microscale* parameter that is
    unexpectedly absent still renders as ``N/A`` and is never silently dropped.

    :param param_specs: Effective ``(attr_name, source, units, fmt)`` specs.
    :return: ``param_specs`` filtered to non-macroscale, non-computed rows.
    """
    macro_names = _macro_param_names()
    return [
        p for p in param_specs
        if p[0] not in macro_names and p[0] not in _COMPUTED
    ]


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

    In the terminal (Rich) output, values that differ from the model defaults
    are highlighted in a distinct color, making non-default runs easy to spot.
    Markdown output (--markdown) is unstyled.

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

    # Formatted default values, for highlighting cells that differ from them.
    # Only the Rich output highlights, so skip this work for Markdown output.
    defaults = None if markdown_out is not None else _default_formatted(
        param_specs, add_names
    )

    if os.path.isfile(path):
        # --- single-file mode ---
        run_code = os.path.splitext(os.path.basename(path))[0]
        data_root = os.path.dirname(path)

        if not no_progress:
            with console.status(f"Loading parameters for {run_code}..."):
                values, has_macro = _load_run_params(
                    data_root, run_code, param_specs, add_names, console
                )
        else:
            values, has_macro = _load_run_params(
                data_root, run_code, param_specs, add_names, console
            )

        if values is None:
            ctx.exit(1)
            return

        # Drop macroscale rows when this file has no macroscale collection.
        effective_specs = param_specs if has_macro else _drop_macro_specs(param_specs)
        nondefault_by_run = (
            None
            if defaults is None
            else {run_code: _nondefault_flags(values, defaults)}
        )
        df = parameters_table(
            {run_code: values},
            effective_specs,
            add_names,
            natural_units,
            [run_code],
            nondefault_by_run=nondefault_by_run,
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
        any_macro = False

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
                vals, has_macro = _load_run_params(
                    path, run_code, param_specs, add_names, console
                )
                if prog is not None:
                    prog.advance(_task)
                if vals is not None:
                    rows[run_code] = vals
                    any_macro = any_macro or has_macro

        if not rows:
            ctx.exit(1)
            return

        # Preserve sort order, skipping any failed runs.  Drop macroscale rows
        # only when no run in the directory has a macroscale collection, so
        # columns still line up when the set is mixed.
        ordered = [rc for rc in run_codes if rc in rows]
        effective_specs = param_specs if any_macro else _drop_macro_specs(param_specs)
        nondefault_by_run = (
            None
            if defaults is None
            else {rc: _nondefault_flags(vals, defaults) for rc, vals in rows.items()}
        )
        df = parameters_table(
            rows,
            effective_specs,
            add_names,
            natural_units,
            ordered,
            nondefault_by_run=nondefault_by_run,
        )

        if markdown_out is not None:
            emit_markdown(
                params_df_to_markdown(df),
                markdown_out,
                console,
            )
        else:
            console.print(params_df_to_rich(df))
