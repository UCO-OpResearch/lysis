"""Experiment initialization from a parameter CSV file.

An :class:`Experiment` is a collection of :class:`~lysis.config.run.Run`
objects that share a common folder.  Each row in the input CSV describes one
Run's parameters; the :meth:`Experiment.from_csv` class method parses the CSV,
resolves any algebraic parameter relationships (via
:mod:`~lysis.config.param_resolver`), validates consistency, creates the
experiment folder on disk, and writes an HDF5 file for every Run.

Typical usage
-------------

**Python API**::

    from lysis.config.experiment import Experiment

    exp = Experiment.from_csv(
        csv_path="my_runs.csv",
        data_root="/data/experiments/",
        name="fiber-radius-sweep",
        description="Effect of fiber radius on lysis time",
    )
    for run in exp.runs:
        print(run.run_code)

**CLI**::

    lysis init-experiment my_runs.csv /data/experiments/ --name fiber-radius-sweep

CSV format
----------

Each **row** describes one Run.  Column headers are Python parameter names
(e.g. ``fiber_radius``, ``pore_size``) — not Fortran names.  Cell values may
include Pint-compatible unit strings (e.g. ``"72.7 nm"``).

Special (non-parameter) columns:

- ``run_code`` — override the auto-generated run code for that row
- ``run_description`` — prose note stored in ``experiment.json``

Any column absent from the CSV is treated as "not provided" and the resolver
will use defaults or solve algebraically from other columns.
"""

import csv
import inspect
import json
import os
from datetime import datetime
from importlib.metadata import version as _pkg_version
from pathlib import Path
from typing import Any

from ..dataio.datastore import DataStore
from .param_resolver import (
    ParameterConflict,
    UnderdeterminedParameters,
    resolve_macro_params,
    resolve_micro_params,
)
from .paramcheck import load_macro_params, load_micro_params
from .parameters import MacroParameters, MicroParameters
from .run import Run


__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


__all__ = ["Experiment"]

# ─── Column classification helpers ───────────────────────────────────────────

# Column names that are Experiment/Run metadata, not parameter names.
_METADATA_COLS = {"run_code", "run_description"}


def _classify_columns(
    headers: list[str],
) -> tuple[set[str], set[str], set[str]]:
    """Partition CSV column names into micro params, macro params, and metadata.

    Uses :func:`inspect.signature` on both parameter classes to determine which
    class each header belongs to.

    :param headers: List of CSV column header strings.
    :returns: ``(micro_names, macro_names, metadata_names)`` — three sets that
        together cover all provided headers.
    :raises ValueError: If any header is not recognised as a micro param, macro
        param, or known metadata column.
    """
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")  # not a CSV column

    # Include dependent (init=False) param names from both classes
    import dataclasses

    micro_dep = {
        f.name
        for f in dataclasses.fields(MicroParameters)
        if f.name not in micro_ind
    }
    macro_dep = {
        f.name
        for f in dataclasses.fields(MacroParameters)
        if f.name not in macro_ind and f.name != "micro_params"
    }

    micro_all = micro_ind | micro_dep
    macro_all = macro_ind | macro_dep

    micro_cols: set[str] = set()
    macro_cols: set[str] = set()
    meta_cols: set[str] = set()
    unknown: list[str] = []

    for h in headers:
        if h in _METADATA_COLS:
            meta_cols.add(h)
        elif h in micro_all:
            micro_cols.add(h)
        elif h in macro_all:
            macro_cols.add(h)
        else:
            unknown.append(h)

    if unknown:
        raise ValueError(
            f"Unrecognised CSV column(s): {unknown}. "
            "Column headers must be Python parameter names from MicroParameters "
            "or MacroParameters, or one of the metadata columns: "
            f"{sorted(_METADATA_COLS)}."
        )

    return micro_cols, macro_cols, meta_cols


# ─── CSV value coercion ───────────────────────────────────────────────────────


def _coerce_csv_value(val_str: str):
    """Convert a CSV string to int, float, bool, or leave as-is for Pint parsing.

    CSV cells that represent pure numbers or booleans arrive as strings from
    :mod:`csv.DictReader`.  Passing them as strings to the parameter classes
    causes ``parse_from_basedict`` to store them as ``str`` instead of the
    expected ``int``/``float``/``bool``, which later breaks h5py attribute
    serialisation.

    Unit-bearing strings (e.g. ``"72.7 nm"``) cannot be parsed as a plain
    number and are returned unchanged so Pint can handle them.

    :param val_str: Raw cell value from the CSV.
    :returns: ``True``/``False`` for boolean strings, ``int`` if the string is
        a plain integer, ``float`` if it is a plain decimal, otherwise the
        original string.
    """
    stripped = val_str.strip()
    if stripped.lower() == "true":
        return True
    if stripped.lower() == "false":
        return False
    try:
        return int(stripped)
    except ValueError:
        pass
    try:
        return float(stripped)
    except ValueError:
        pass
    return stripped


# ─── Run-code generation ──────────────────────────────────────────────────────


def _make_run_code(base_ts: str, row_index: int) -> str:
    """Generate a unique run code for *row_index* within a single CSV parse.

    Appends a zero-padded two-digit row index to *base_ts* (the timestamp) to
    prevent collisions when all rows are processed within the same minute.

    :param base_ts: Timestamp string from ``datetime.now().strftime("%Y-%m-%d-%H%M")``.
    :param row_index: Zero-based row index within the CSV.
    :returns: Run code string, e.g. ``"2026-04-07-1422-00"``.
    """
    return f"{base_ts}-{row_index:02d}"


# ─── Experiment class ─────────────────────────────────────────────────────────


class Experiment:
    """A collection of :class:`~lysis.config.run.Run` objects in a shared folder.

    Normally created via :meth:`from_csv` rather than directly.

    :param name: Human-readable experiment identifier; used as the folder name.
    :param data_root: Parent directory under which the experiment folder is
        created.
    :param description: Optional prose description stored in
        ``experiment.json``.
    """

    def __init__(
        self,
        name: str,
        data_root: str | os.PathLike,
        description: str = "",
    ):
        self._name = name
        self._data_root = str(data_root)
        self._description = description
        self._runs: list[Run] = []
        self._created: str = datetime.now().isoformat(timespec="seconds")

    # ─── Properties ────────────────────────────────────────────────────────

    @property
    def name(self) -> str:
        """The experiment name (also the folder name)."""
        return self._name

    @property
    def path(self) -> str:
        """Absolute path to the experiment folder."""
        return os.path.join(self._data_root, self._name)

    @property
    def runs(self) -> list[Run]:
        """List of :class:`~lysis.config.run.Run` objects, one per CSV row."""
        return list(self._runs)

    @property
    def description(self) -> str:
        """Prose description of the experiment."""
        return self._description

    # ─── Serialisation ─────────────────────────────────────────────────────

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable representation of this Experiment.

        Suitable for writing ``experiment.json``.

        :returns: Dict with keys ``name``, ``description``, ``created``,
            ``lysis_version``, and ``runs`` (list of run summaries).
        :rtype: dict[str, Any]
        """
        try:
            lysis_ver = _pkg_version("lysis")
        except Exception:
            lysis_ver = "unknown"

        return {
            "name": self._name,
            "description": self._description,
            "created": self._created,
            "lysis_version": lysis_ver,
            "runs": [
                {
                    "run_code": run.run_code,
                    "row_index": i,
                    "description": getattr(run, "_description", ""),
                }
                for i, run in enumerate(self._runs)
            ],
        }

    # ─── Class method constructor ─────────────────────────────────────────

    @classmethod
    def from_csv(
        cls,
        csv_path: str | os.PathLike,
        data_root: str | os.PathLike,
        name: str | None = None,
        description: str = "",
        dry_run: bool = False,
    ) -> "Experiment":
        """Parse a parameter CSV and initialise an Experiment.

        Each row in the CSV describes one Run.  The method:

        1. Parses the CSV with :mod:`csv.DictReader`.
        2. Classifies columns as micro params, macro params, or metadata.
        3. For each row: resolves algebraic parameter relationships, validates
           consistency, and constructs :class:`~.parameters.MicroParameters` /
           :class:`~.parameters.MacroParameters` objects.
        4. Collects **all** row errors before raising, so researchers see every
           problem at once.
        5. If ``dry_run`` is ``False``: creates the experiment folder,
           writes ``experiment.json``, and calls
           :meth:`~lysis.dataio.datastore.DataStore.create` +
           :meth:`~lysis.dataio.datastore.DataStore.initialize_macroscale`
           for each Run.

        :param csv_path: Path to the parameter CSV file.
        :type csv_path: str | os.PathLike
        :param data_root: Parent directory for the new experiment folder.
        :type data_root: str | os.PathLike
        :param name: Experiment name (folder name).  Defaults to the CSV
            filename stem.
        :type name: str | None
        :param description: Optional prose description stored in
            ``experiment.json``.
        :type description: str
        :param dry_run: If ``True``, validate and resolve parameters but do
            **not** create any files or folders.
        :type dry_run: bool
        :returns: :class:`Experiment` instance with all Runs populated.
        :rtype: Experiment
        :raises FileExistsError: If the experiment folder already exists (only
            when ``dry_run=False``).
        :raises ValueError: If any CSV row contains unrecognised column names.
        :raises ParameterConflict: If any CSV row contains mutually inconsistent
            parameter values.  All row errors are collected before raising.
        :raises ValueError: If ``load_micro_params`` or ``load_macro_params``
            finds missing required parameters or dependent-value mismatches.
        """
        csv_path = Path(csv_path)
        if name is None:
            name = csv_path.stem

        exp = cls(name, data_root, description)

        # ── 1. Parse CSV ────────────────────────────────────────────────
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file is empty: {csv_path}")
            micro_cols, macro_cols, meta_cols = _classify_columns(
                list(reader.fieldnames)
            )
            rows = list(reader)

        if not rows:
            raise ValueError(f"CSV file has no data rows: {csv_path}")

        # ── 2. Resolve parameters for every row (collect all errors) ────
        base_ts = datetime.now().strftime("%Y-%m-%d-%H%M")
        resolved_rows: list[dict] = []  # {micro_params, macro_params, run_code, run_description}
        errors: list[str] = []

        for row_idx, row in enumerate(rows):
            micro_provided = {
                k: _coerce_csv_value(v)
                for k, v in row.items()
                if k in micro_cols and v != ""
            }
            macro_provided = {
                k: _coerce_csv_value(v)
                for k, v in row.items()
                if k in macro_cols and v != ""
            }
            run_code = row.get("run_code", "").strip() or _make_run_code(base_ts, row_idx)
            run_desc = row.get("run_description", "").strip()

            try:
                resolved_micro = resolve_micro_params(micro_provided)
                micro_params = load_micro_params(resolved_micro)
            except (ParameterConflict, ValueError) as exc:
                errors.append(f"Row {row_idx + 1} (run_code={run_code!r}) micro: {exc}")
                continue

            try:
                resolved_macro = resolve_macro_params(
                    macro_provided, micro_params.to_basedict()
                )
                macro_params = load_macro_params(resolved_macro, micro_params)
            except (ParameterConflict, ValueError) as exc:
                errors.append(f"Row {row_idx + 1} (run_code={run_code!r}) macro: {exc}")
                continue

            resolved_rows.append(
                {
                    "micro_params": micro_params,
                    "macro_params": macro_params,
                    "run_code": run_code,
                    "run_description": run_desc,
                }
            )

        if errors:
            raise ValueError(
                f"Parameter errors in {csv_path.name} "
                f"({len(errors)} row(s)):\n"
                + "\n".join(errors)
            )

        # ── 3. dry_run stops here ────────────────────────────────────────
        if dry_run:
            # Populate _runs with lightweight Run objects (no HDF5)
            for rr in resolved_rows:
                run = Run(exp._data_root, rr["run_code"])
                run.micro_params = rr["micro_params"]
                run.macro_params = rr["macro_params"]
                run._description = rr["run_description"]
                exp._runs.append(run)
            return exp

        # ── 4. Create experiment folder ──────────────────────────────────
        exp_path = exp.path
        os.makedirs(exp_path, exist_ok=False)

        # ── 5. Write experiment.json ─────────────────────────────────────
        # (populate _runs first so to_dict() includes them)
        for rr in resolved_rows:
            run = Run(exp_path, rr["run_code"])
            run.micro_params = rr["micro_params"]
            run.macro_params = rr["macro_params"]
            run._description = rr["run_description"]
            exp._runs.append(run)

        json_path = os.path.join(exp_path, "experiment.json")
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(exp.to_dict(), fh, indent=2)

        # ── 6. Create HDF5 files for each Run ───────────────────────────
        for run in exp._runs:
            ds = DataStore.create(run.run_code, exp_path, run.micro_params)
            ds.initialize_macroscale(run.macro_params)
            ds.close()

        return exp
