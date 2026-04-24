"""Experiment initialization from a parameter CSV file.

An :class:`Experiment` is a collection of :class:`~lysis.config.run.Run`
objects that share a common folder.  Each column in the input CSV describes one
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

Each **column** describes one Run.  The first column contains parameter names
(e.g. ``fiber_radius``, ``pore_size``) — not Fortran names.  Each additional
column header is the run code for that Run; remaining cells are the
corresponding parameter values, which may include Pint-compatible unit strings
(e.g. ``"72.7 nm"``).

Special (non-parameter) rows:

- ``run_description`` — prose note stored in ``experiment.json`` (one cell per run)

The column header (row 0) is the ``run_code``; if left blank an auto-generated
timestamp-based code is used.  Any row absent from the CSV is treated as
"not provided" and the resolver will use defaults or solve algebraically.
"""

import csv
import inspect
import json
import math
import os
from datetime import datetime
from importlib.metadata import version as _pkg_version
from pathlib import Path
from typing import Any

from ..dataio.datastore import DataStore
from ..dataio.fileops import _params_json_default
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
    """Partition parameter names into micro params, macro params, and metadata.

    Uses :func:`inspect.signature` on both parameter classes to determine which
    class each name belongs to.  In the transposed CSV format these names come
    from the first column of data rows (not from column headers).

    :param headers: List of parameter name strings (first column of data rows).
    :returns: ``(micro_names, macro_names, metadata_names)`` — three sets that
        together cover all provided names.
    :raises ValueError: If any name is not recognised as a micro param, macro
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
        # Populated by ``Experiment.load``; ``None`` for Experiments created
        # via :meth:`from_csv` (where ``_name`` and the folder name are
        # guaranteed to match by construction).
        self._source_path: str | None = None

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

        def _serialisable_macro(macro_params):
            """Return macro_params as a JSON-safe dict (NaN values excluded)."""
            if macro_params is None:
                return None
            return {
                k: v
                for k, v in macro_params.to_basedict().items()
                if not (isinstance(v, float) and math.isnan(v))
            }

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
                    "macro_params": _serialisable_macro(run.macro_params),
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

        Each column in the CSV describes one Run.  The first column contains
        parameter names; each subsequent column header is the run code and
        the cells below it are values for that Run.  The method:

        1. Parses the transposed CSV with :mod:`csv.reader`.
        2. Classifies parameter names (first column) as micro params, macro
           params, or metadata.
        3. For each run column: resolves algebraic parameter relationships, validates
           consistency, and constructs :class:`~.parameters.MicroParameters` /
           :class:`~.parameters.MacroParameters` objects.
        4. Collects **all** row errors before raising, so researchers see every
           problem at once.
        5. If ``dry_run`` is ``False``: creates the experiment folder,
           writes ``experiment.json``, and calls
           :meth:`~lysis.dataio.datastore.DataStore.create` for each Run,
           producing an HDF5 file with microscale parameters and empty
           microscale dataset stubs.

        .. note::
            Macroscale parameters are stored in ``experiment.json`` and in
            each :class:`~lysis.config.run.Run` object but are **not** written
            to the HDF5 file at this stage.
            :meth:`~lysis.dataio.datastore.DataStore.initialize_macroscale`
            must be called after all microscale Simulations have completed;
            it reads the microscale results to compute ``forced_unbind``
            and then writes the full macroscale structure to the HDF5.

        :param csv_path: Path to the parameter CSV file (transposed format:
            first column = parameter names, each additional column = one Run).
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
        :raises ValueError: If any parameter name in the first column is not recognised.
        :raises ParameterConflict: If any run column contains mutually inconsistent
            parameter values.  All run errors are collected before raising.
        :raises ValueError: If ``load_micro_params`` or ``load_macro_params``
            finds missing required parameters or dependent-value mismatches.
        """
        csv_path = Path(csv_path)
        if name is None:
            name = csv_path.stem

        exp = cls(name, data_root, description)

        # ── 1. Parse CSV (transposed: rows = params, columns = runs) ────
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            all_rows = list(reader)

        if not all_rows:
            raise ValueError(f"CSV file is empty: {csv_path}")

        # Row 0: ["parameter", run_code_0, run_code_1, ...]
        header_row = all_rows[0]
        n_runs = len(header_row) - 1
        if n_runs <= 0:
            raise ValueError(f"CSV has no run columns: {csv_path}")

        raw_run_codes = [cell.strip() for cell in header_row[1:]]

        # Data rows: [param_name, value_for_run_0, value_for_run_1, ...]
        data_rows = [r for r in all_rows[1:] if any(cell.strip() for cell in r)]
        if not data_rows:
            raise ValueError(f"CSV file has no data rows: {csv_path}")

        # Pad rows that are shorter than expected
        for r in data_rows:
            while len(r) <= n_runs:
                r.append("")

        # Skip rows with a blank parameter name
        data_rows = [r for r in data_rows if r[0].strip()]

        param_names = [r[0].strip() for r in data_rows]
        micro_cols, macro_cols, meta_cols = _classify_columns(param_names)

        # Build one dict per run column
        rows = [
            {r[0].strip(): r[col_idx + 1] for r in data_rows}
            for col_idx in range(n_runs)
        ]

        # ── 2. Resolve parameters for every row (collect all errors) ────
        base_ts = datetime.now().strftime("%Y-%m-%d-%H%M")
        resolved_rows: list[dict] = []  # {micro_params, macro_params, run_code, run_description}
        errors: list[str] = []

        for col_idx, row in enumerate(rows):
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
            run_code = raw_run_codes[col_idx] or _make_run_code(base_ts, col_idx)
            run_desc = row.get("run_description", "").strip()

            try:
                resolved_micro = resolve_micro_params(micro_provided)
                micro_params = load_micro_params(resolved_micro)
            except (ParameterConflict, ValueError) as exc:
                errors.append(f"Run {col_idx + 1} (run_code={run_code!r}) micro: {exc}")
                continue

            try:
                resolved_macro = resolve_macro_params(
                    macro_provided, micro_params.to_basedict()
                )
                macro_params = load_macro_params(resolved_macro, micro_params)
            except (ParameterConflict, ValueError) as exc:
                errors.append(f"Run {col_idx + 1} (run_code={run_code!r}) macro: {exc}")
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
                f"({len(errors)} run(s)):\n"
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
            json.dump(
                exp.to_dict(), fh, indent=2, default=_params_json_default
            )

        # ── 6. Create HDF5 files for each Run ───────────────────────────
        # Only the microscale structure is written here. initialize_macroscale()
        # must be called separately after microscale Simulations complete,
        # because it requires microscale output to compute forced_unbind.
        for run in exp._runs:
            ds = DataStore.create(run.run_code, exp_path, run.micro_params)
            ds.close()

        return exp

    # ─── Renaming ──────────────────────────────────────────────────────────

    def _write_json(self) -> None:
        """Write the current state of this Experiment to ``experiment.json``."""
        json_path = os.path.join(self.path, "experiment.json")
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2, default=_params_json_default)

    def rename(self, new_name: str) -> str:
        """Rename this Experiment's folder and update ``experiment.json``.

        Renames the on-disk folder from ``{data_root}/{name}`` to
        ``{data_root}/{new_name}``, updates ``self._name`` and the
        ``os_path`` of every contained :class:`~lysis.config.run.Run`, and
        rewrites ``experiment.json`` with the new name.

        If the Experiment was loaded from disk via :meth:`load` and the
        folder basename no longer matches ``experiment.json``'s ``name``
        field, this method refuses to rename rather than silently
        widening the divergence.  Fix either the folder name or the
        ``name`` field manually before retrying.

        :param new_name: The new experiment name (becomes the new folder
            name under ``data_root``).
        :type new_name: str
        :return: The previous experiment name.
        :rtype: str
        :raises ValueError: If ``new_name`` equals the current name, or
            if the loaded folder basename does not match
            ``experiment.json``'s ``name`` field.
        :raises FileNotFoundError: If the current experiment folder does
            not exist.
        :raises FileExistsError: If a folder already exists at the
            destination path.
        """
        if new_name == self._name:
            raise ValueError(
                f"new_name is identical to current name: {new_name!r}"
            )

        if self._source_path is not None:
            loaded_name = os.path.basename(self._source_path)
            if loaded_name != self._name:
                raise ValueError(
                    f"Folder name ({loaded_name!r}) and experiment.json "
                    f"name ({self._name!r}) disagree. Refusing to rename; "
                    f"fix the folder name or the 'name' field in "
                    f"experiment.json before retrying."
                )

        old_name = self._name
        old_path = self.path
        new_path = os.path.join(self._data_root, new_name)

        if not os.path.isdir(old_path):
            raise FileNotFoundError(f"Experiment folder not found: {old_path}")
        if os.path.exists(new_path):
            raise FileExistsError(f"Path already exists: {new_path}")

        os.rename(old_path, new_path)

        self._name = new_name
        self._source_path = new_path
        for run in self._runs:
            run.os_path = new_path

        self._write_json()
        return old_name

    def rename_run(self, old_run_code: str, new_run_code: str) -> None:
        """Rename a Run belonging to this Experiment.

        Locates the :class:`~lysis.config.run.Run` with ``run_code ==
        old_run_code``, validates that neither the experiment manifest
        nor the on-disk HDF5 files already use ``new_run_code``, calls
        :meth:`~lysis.config.run.Run.rename` on the target (which renames
        the HDF5 file and records the history in the ``renamed_from``
        attribute), verifies the rename completed, then rewrites
        ``experiment.json``.

        The json rewrite is deliberately the *second* step and it
        confirms the HDF5 rename actually landed on disk before touching
        the manifest, so a failure between the two steps surfaces as an
        error rather than a silently inconsistent experiment.

        :param old_run_code: The current ``run_code`` of the Run.
        :type old_run_code: str
        :param new_run_code: The new ``run_code``.
        :type new_run_code: str
        :raises KeyError: If no Run in this Experiment has
            ``run_code == old_run_code``.
        :raises ValueError: If ``new_run_code`` equals ``old_run_code``,
            or if another Run in this Experiment already uses
            ``new_run_code``.
        :raises FileNotFoundError: If the source HDF5 file does not exist.
        :raises FileExistsError: If an HDF5 file already exists at the
            destination path.
        :raises RuntimeError: If :meth:`Run.rename` returned but the
            HDF5 file is not at its expected new location.
        """
        if new_run_code == old_run_code:
            raise ValueError(
                f"new_run_code is identical to old_run_code: {new_run_code!r}"
            )

        target = None
        for run in self._runs:
            if run.run_code == old_run_code:
                target = run
                break
        if target is None:
            raise KeyError(
                f"No run with run_code={old_run_code!r} in experiment {self._name!r}"
            )
        if any(r.run_code == new_run_code for r in self._runs):
            raise ValueError(
                f"run_code={new_run_code!r} is already in use in experiment "
                f"{self._name!r}"
            )

        old_h5 = os.path.join(target.os_path, f"{old_run_code}.h5")
        new_h5 = os.path.join(target.os_path, f"{new_run_code}.h5")
        if not os.path.isfile(old_h5):
            raise FileNotFoundError(f"HDF5 file not found: {old_h5}")
        if os.path.exists(new_h5):
            raise FileExistsError(f"HDF5 file already exists: {new_h5}")

        # Step 1: rename the HDF5 file (Run.rename re-checks file state).
        target.rename(new_run_code)

        # Step 2 verifies Step 1 landed before we touch experiment.json.
        if not os.path.isfile(new_h5) or os.path.exists(old_h5):
            raise RuntimeError(
                f"Run.rename reported success but the HDF5 file is not "
                f"where expected (missing: {new_h5}, lingering: {old_h5}). "
                f"experiment.json was not updated."
            )

        self._write_json()

    @classmethod
    def load(cls, experiment_path: str | os.PathLike) -> "Experiment":
        """Load an existing Experiment from its folder on disk.

        Reads ``experiment.json`` to reconstruct all :class:`~lysis.config.run.Run`
        objects, including their macro- and microscale parameters.  Microscale
        parameters are read from each run's HDF5 file; macroscale parameters are
        read from ``experiment.json`` (where they were stored by :meth:`from_csv`).

        :param experiment_path: Path to the experiment folder (the directory
            that contains ``experiment.json`` and the run HDF5 files).
        :type experiment_path: str | os.PathLike
        :returns: :class:`Experiment` instance with all :class:`~lysis.config.run.Run`
            objects populated.
        :rtype: Experiment
        :raises FileNotFoundError: If ``experiment.json`` is not found.
        :raises KeyError: If ``experiment.json`` is missing required fields.
        :raises ValueError: If macro_params stored in ``experiment.json`` are
            inconsistent with the micro_params in the corresponding HDF5 file.
        """
        experiment_path = Path(experiment_path).resolve()
        json_path = experiment_path / "experiment.json"
        if not json_path.exists():
            raise FileNotFoundError(
                f"experiment.json not found in {experiment_path}. "
                "Is this a valid experiment folder?"
            )

        with open(json_path, encoding="utf-8") as fh:
            meta = json.load(fh)

        name = meta["name"]
        description = meta.get("description", "")
        exp = cls(name, str(experiment_path.parent), description)
        exp._created = meta.get("created", exp._created)
        exp._source_path = str(experiment_path)

        for run_meta in meta["runs"]:
            run_code = run_meta["run_code"]
            macro_params_dict = run_meta.get("macro_params") or {}
            run_desc = run_meta.get("description", "")

            # Micro parameters come from the HDF5 file
            ds = DataStore(run_code, str(experiment_path), mode="r")
            micro_params = ds.micro_params
            ds.close()

            # Macro parameters come from experiment.json; forced_unbind is
            # absent (NaN was excluded when writing) so load_macro_params
            # will fill in the nan default.  initialize_macroscale() will
            # compute the real value when called.
            macro_params = load_macro_params(macro_params_dict, micro_params)

            run = Run(str(experiment_path), run_code)
            run.micro_params = micro_params
            run.macro_params = macro_params
            run._description = run_desc
            exp._runs.append(run)

        return exp
