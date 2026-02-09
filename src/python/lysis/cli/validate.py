"""``lysis validate`` — validate simulation data against a specification."""

import sys

import click
import numpy as np

from lysis.cli import cli
from lysis.cli.convert import _resolve_spec
from lysis.data.dataspec import dataspec, check_dataset_spec


@cli.command()
@click.argument("data_path", type=click.Path(exists=True))
@click.option(
    "-s",
    "--spec",
    required=True,
    help="Spec version or tag to validate against (e.g. hdf5, v2.0.0).",
)
@click.option(
    "-c",
    "--collections",
    default=None,
    help="Comma-separated collection names (default: all in spec).",
)
@click.option(
    "--file-code",
    default="",
    help="File code for Fortran file naming.",
)
@click.pass_context
def validate(ctx, data_path, spec, collections, file_code):
    """Validate simulation data against a specification.

    Reads data from DATA_PATH and checks each dataset against the expected
    dtype and shape defined in the specification.  Returns exit code 0 if
    all datasets pass, 1 if any fail.

    \b
    Examples:
        lysis validate output.h5 -s hdf5
        lysis validate ./fortran_data/ -s fortran --file-code "_run01.dat"
        lysis validate output.h5 -s v2.0.0 -c microscale_out
    """
    from lysis.data.fileops import read_data_collection

    console = ctx.obj["console"]
    verbose = ctx.obj["verbose"]

    # Resolve spec version
    try:
        spec_version = _resolve_spec(spec)
    except click.BadParameter as e:
        console.print(f"[red]Error:[/red] {e.format_message()}")
        ctx.exit(1)
        return

    # Determine collections
    if collections:
        collection_names = [c.strip() for c in collections.split(",")]
        for name in collection_names:
            if name not in dataspec[spec_version]:
                available = list(dataspec[spec_version].keys())
                console.print(
                    f"[red]Error:[/red] Unknown collection '{name}'. "
                    f"Available: {', '.join(available)}"
                )
                ctx.exit(1)
                return
    else:
        collection_names = list(dataspec[spec_version].keys())

    spec_collections = [dataspec[spec_version][n] for n in collection_names]
    file_codes = [file_code] * len(collection_names)

    console.print(f"Validating against [bold]{spec_version}[/bold]...\n")

    # Read data
    try:
        data = read_data_collection(data_path, spec_collections, file_codes)
    except FileNotFoundError as e:
        console.print(f"[red]Error:[/red] Data file not found: {e}")
        ctx.exit(1)
        return
    except KeyError as e:
        console.print(f"[red]Error:[/red] Missing data: {e}")
        ctx.exit(1)
        return

    params = data.get("params")

    # Validate each dataset
    passed = 0
    failed = 0
    total = 0

    for name in collection_names:
        collection = dataspec[spec_version][name]
        console.print(f"[bold]{name}:[/bold]")

        for dataset_name, dataset_spec in collection.data.items():
            total += 1

            if dataset_name not in data:
                console.print(
                    f"  {dataset_name:<30} [red]MISSING[/red]"
                )
                failed += 1
                continue

            dataset_data = data[dataset_name]

            # Handle per-simulation lists
            if isinstance(dataset_data, list):
                if not dataset_data:
                    console.print(
                        f"  {dataset_name:<30} [red]FAIL[/red]  empty list"
                    )
                    failed += 1
                    continue

                all_ok = all(
                    check_dataset_spec(d, dataset_spec, params=params)
                    for d in dataset_data
                )
                sample = dataset_data[0]
                count = len(dataset_data)
                suffix = f"  [{count} sim{'s' if count != 1 else ''}]"
            else:
                all_ok = check_dataset_spec(
                    dataset_data, dataset_spec, params=params
                )
                sample = dataset_data
                suffix = ""

            dtype_str = str(sample.dtype)
            shape_str = str(sample.shape)

            if all_ok:
                console.print(
                    f"  {dataset_name:<30} [green]OK[/green]    "
                    f"{dtype_str:<10} {shape_str}{suffix}"
                )
                passed += 1
            else:
                # Build failure detail
                detail_parts = []
                expected_dtype = np.dtype(dataset_spec.dtype)
                if not np.can_cast(sample.dtype, expected_dtype, casting="same_kind"):
                    detail_parts.append(
                        f"dtype: got {sample.dtype}, expected {expected_dtype}"
                    )
                else:
                    detail_parts.append(f"shape mismatch: {shape_str}")

                detail = "; ".join(detail_parts)
                console.print(
                    f"  {dataset_name:<30} [red]FAIL[/red]  {detail}{suffix}"
                )
                failed += 1

        console.print()  # Blank line between collections

    # Summary
    if failed == 0:
        console.print(
            f"[green]Result: {passed}/{total} datasets passed[/green]"
        )
    else:
        console.print(
            f"[red]Result: {passed}/{total} passed, {failed} failed[/red]"
        )
        ctx.exit(1)
