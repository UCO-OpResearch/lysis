"""``lysis convert`` — convert simulation data between specification formats."""

import json
import sys

import click

from lysis.cli import cli
from lysis.dataio.dataspec import dataspec, tags


def _resolve_spec(spec_str):
    """Resolve a spec name or tag to a version string.

    :param spec_str: Specification version or tag alias
        (e.g. ``"fortran"``, ``"hdf5"``, ``"v1.99.0"``)
    :type spec_str: str
    :return: Resolved version string (e.g. ``"v1.99.0"``)
    :rtype: str
    :raises click.BadParameter: If the spec is not recognized
    """
    resolved = spec_str
    while resolved in tags:
        resolved = tags[resolved]
    if resolved not in dataspec:
        available = sorted(k for k in dataspec if not k in tags and k != resolved)
        tag_list = [f"{k} -> {v}" for k, v in tags.items()]
        raise click.BadParameter(
            f"Unknown spec '{spec_str}'.\n"
            f"  Versions: {', '.join(available)}\n"
            f"  Tags: {', '.join(tag_list)}"
        )
    return resolved


@cli.command()
@click.argument("input_path", type=click.Path(exists=True))
@click.argument("output_path", type=click.Path())
@click.option(
    "-f",
    "--from",
    "input_spec",
    required=True,
    help="Input spec version or tag (e.g. fortran, hdf5, v1.99.0, v2.0.0).",
)
@click.option(
    "-t",
    "--to",
    "output_spec",
    required=True,
    help="Output spec version or tag.",
)
@click.option(
    "-c",
    "--collections",
    default=None,
    help="Comma-separated collection names (default: all in input spec).",
)
@click.option(
    "--file-code",
    default="",
    help=(
        "File code(s) for Fortran file naming (e.g. '_PLG2_tPA01_TB-xiii.dat')."
        "Can be blank (default), a single string that will be used for all collections, "
        "or a comma-separated list matching the list of collections."
        "Example: --file-code _PLG2_tPA01_TB-xiii,_PLG2_tPA01_TB-xiii,_TB-xiii__21_105"
    ),
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be converted without writing output.",
)
@click.option(
    "--no-progress",
    is_flag=True,
    default=False,
    help="Suppress progress indicators.",
)
@click.option(
    "--param-override",
    "param_overrides",
    multiple=True,
    metavar="KEY=VALUE",
    help=(
        "Override a parameter value during validation (repeatable). "
        "KEY must be a member of MicroParameters or MacroParameters. "
        "VALUE can be a scalar or a string that will parse as a Pint Quantity. "
        "Example: --param-override total_molecules=100 --param-override pore_size=1.0135um"
    ),
)
@click.option(
    "--param-alias",
    "param_aliases",
    multiple=True,
    metavar="KEY=FORTRAN_NAME",
    help=(
        "Alias a Fortran variable name to a Python parameter (repeatable). "
        "KEY is a Python parameter name; FORTRAN_NAME is the name in the log file. "
        "Example: --param-alias bind_rate_tPA=kon"
    ),
)
@click.pass_context
def convert(
    ctx,
    input_path,
    output_path,
    input_spec,
    output_spec,
    collections,
    file_code,
    dry_run,
    no_progress,
    param_overrides,
    param_aliases,
):
    """Convert simulation data between specification formats.

    Reads data from INPUT_PATH in the --from format, converts it, and writes
    it to OUTPUT_PATH in the --to format.

    \b
    Examples:
        lysis convert ./fortran_data/ output.h5 -f fortran -t hdf5
        lysis convert input.h5 ./out/ -f hdf5 -t fortran --file-code "_run01.dat"
        lysis convert input.h5 output.h5 -f v1.99.0 -t v2.0.0 --dry-run
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --param-override total_molecules=100
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --param-override pore_size=1.0135um
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --param-alias bind_rate_tPA=kon
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --micro-log micro_PLG2_tPA01.txt
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --macro-log macro_TB-xiii__21_105.txt
        lysis convert ./data/ out.h5 -f fortran -t hdf5 --micro-log micro.txt --macro-log macro.txt
    """
    from lysis.dataio.dataconvert import convert_data
    from lysis.dataio.fileops import read_data_collection, write_data_collection

    console = ctx.obj["console"]
    verbose = ctx.obj["verbose"]

    # Resolve spec versions
    try:
        in_version = _resolve_spec(input_spec)
        out_version = _resolve_spec(output_spec)
    except click.BadParameter as e:
        console.print(f"[red]Error:[/red] {e.format_message()}")
        ctx.exit(1)
        return

    # Determine collections to convert
    if collections:
        collection_names = [c.strip() for c in collections.split(",")]
        for name in collection_names:
            if name not in dataspec[in_version]:
                available = list(dataspec[in_version].keys())
                console.print(
                    f"[red]Error:[/red] Unknown collection '{name}'. "
                    f"Available: {', '.join(available)}"
                )
                ctx.exit(1)
                return
    else:
        collection_names = list(dataspec[in_version].keys())

    in_collections = [dataspec[in_version][n] for n in collection_names]
    out_collections = [dataspec[out_version][n] for n in collection_names]
    if file_code:
        file_codes = [f.strip() for f in file_code.split(",")]
        if len(file_codes) == 1:
            file_codes = file_codes * len(collection_names)
        elif not len(file_codes) == len(collection_names):
            console.print(
                f"[red]Error:[/red] --file-code must be none, one, or one-per-collection, got: {item!r}"
            )
            ctx.exit(1)
            return
    else:
        file_codes = [""] * len(collection_names)

    # Validate that a conversion path exists before reading any data
    if in_version != out_version:
        from lysis.dataio.dataconvert import conversion_paths

        if (in_version, out_version) not in conversion_paths:
            console.print(
                f"[red]Error:[/red] No conversion path from "
                f"{in_version} to {out_version}."
            )
            ctx.exit(1)
            return
        path = conversion_paths[in_version, out_version]
    else:
        path = [in_version]

    if verbose or dry_run:
        console.print(f"Converting: {in_version} -> {out_version}")
        if len(path) > 2:
            console.print(f"Path: {' -> '.join(path)}")
        console.print(f"Collections: {', '.join(collection_names)}")
        console.print(f"Input:  {input_path}")
        console.print(f"Output: {output_path}")
        if dry_run:
            console.print("[yellow]Dry run — no files will be written.[/yellow]")

    # Parse --param-override KEY=VALUE pairs
    overrides = None
    if param_overrides:
        overrides = {}
        for item in param_overrides:
            if "=" not in item:
                console.print(
                    f"[red]Error:[/red] --param-override must be KEY=VALUE, got: {item!r}"
                )
                ctx.exit(1)
                return
            key, raw = item.split("=", 1)
            try:
                overrides[key] = json.loads(raw)
            except json.JSONDecodeError:
                overrides[key] = raw  # treat as bare string

    # Parse --param-alias KEY=FORTRAN_NAME pairs
    aliases = None
    if param_aliases:
        aliases = {}
        for item in param_aliases:
            if "=" not in item:
                console.print(
                    f"[red]Error:[/red] --param-alias must be KEY=FORTRAN_NAME, got: {item!r}"
                )
                ctx.exit(1)
                return
            key, fortran_name = item.split("=", 1)
            aliases[key] = fortran_name

    from contextlib import nullcontext

    from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

    _pctx = (
        Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            TimeElapsedColumn(),
            console=console,
        )
        if not no_progress
        else nullcontext()
    )

    with _pctx as prog:
        _task = (
            prog.add_task("Reading data...", total=None) if prog is not None else None
        )

        # Read
        try:
            data = read_data_collection(input_path, in_collections, file_codes)
        except FileNotFoundError as e:
            console.print(f"[red]Error:[/red] Input file not found: {e}")
            ctx.exit(1)
            return
        except KeyError as e:
            console.print(f"[red]Error:[/red] Missing data in input: {e}")
            ctx.exit(1)
            return

        # Apply param_overrides: merge into all params sub-dicts
        if overrides:
            for key, value in overrides.items():
                for section in data["params"].values():
                    if isinstance(section, dict):
                        section[key] = value

        # Apply param_aliases: rename keys in all params sub-dicts
        if aliases:
            for py_name, fort_name in aliases.items():
                fort_lower = fort_name.lower()
                for section in data["params"].values():
                    if isinstance(section, dict) and fort_lower in section:
                        section[py_name] = section.pop(fort_lower)

        if verbose:
            n_datasets = sum(1 for k in data if k != "params")
            console.print(f"Read {n_datasets} dataset(s) from input.")

        if prog is not None:
            prog.update(_task, description="Converting data...")

        # Convert
        try:
            converted = convert_data(data, in_version, out_version)
        except ValueError as e:
            console.print(f"[red]Error:[/red] Conversion error: {e}")
            ctx.exit(1)
            return
        except NotImplementedError as e:
            console.print(f"[red]Error:[/red] Conversion not implemented: {e}")
            ctx.exit(1)
            return
        except (OverflowError, TypeError) as e:
            console.print(f"[red]Error:[/red] Type conversion failed: {e}")
            ctx.exit(1)
            return

        # Write
        if dry_run:
            console.print("[green]Dry run complete. No files written.[/green]")
            return

        if prog is not None:
            prog.update(_task, description="Writing output...")

        try:
            write_data_collection(converted, output_path, out_collections, file_codes)
        except (OSError, TypeError) as e:
            console.print(f"[red]Error:[/red] Failed to write output: {e}")
            ctx.exit(1)
            return

        if prog is not None:
            prog.update(_task, description="[green]Done.[/green]")

    console.print(
        f"[green]Converted {in_version} -> {out_version} "
        f"({', '.join(collection_names)})[/green]"
    )
