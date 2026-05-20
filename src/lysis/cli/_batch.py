"""Shared helpers for batch (experiment-folder) execution in CLI commands.

Used by ``lysis run-micro`` and ``lysis run-macro`` when a directory is passed
instead of a single HDF5 file.
"""

import json
from pathlib import Path

import click


__all__ = ["resolve_hdf5_paths"]


def resolve_hdf5_paths(directory: "str | Path") -> "list[Path]":
    """Return the ordered list of HDF5 :class:`~pathlib.Path` objects in *directory*.

    Resolution strategy:

    1. If ``experiment.json`` exists in *directory*, the run codes are read
       from the JSON (in CSV row order) and each corresponding ``{run_code}.h5``
       path is returned.  Any missing HDF5 file raises
       :exc:`FileNotFoundError` immediately.
    2. Otherwise, glob for ``*.h5`` files and return them sorted by name.

    :param directory: Path to a folder containing HDF5 simulation files.
    :type directory: str or Path
    :returns: Non-empty list of HDF5 :class:`~pathlib.Path` objects.
    :rtype: list[Path]
    :raises click.UsageError: If the directory contains no ``.h5`` files
        (glob mode) or if the directory argument is not actually a directory.
    :raises FileNotFoundError: If a run listed in ``experiment.json`` has no
        corresponding HDF5 file on disk.
    """
    directory = Path(directory)
    json_path = directory / "experiment.json"

    if json_path.exists():
        with open(json_path, encoding="utf-8") as fh:
            meta = json.load(fh)
        paths: list[Path] = []
        for run_meta in meta.get("runs", []):
            run_code = run_meta["run_code"]
            p = directory / f"{run_code}.h5"
            if not p.exists():
                raise FileNotFoundError(
                    f"Expected HDF5 file not found: {p}"
                )
            paths.append(p)
        if not paths:
            raise click.UsageError(
                f"experiment.json in {directory} lists no runs."
            )
        return paths

    # No experiment.json — fall back to sorted glob
    paths = sorted(directory.glob("*.h5"))
    if not paths:
        raise click.UsageError(f"No .h5 files found in {directory}")
    return paths
