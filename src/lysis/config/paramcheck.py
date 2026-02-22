"""Parameter validation and Fortran log verification.

Provides strict parameter loading and Fortran log parsing/verification on top
of :meth:`~.parameters.Parameters.parse_from_basedict`.  All functions here
are **optional** — existing code that calls ``parse_from_basedict()`` directly
is unchanged.

Strict loading
--------------

:func:`load_micro_params` and :func:`load_macro_params` raise
:exc:`ValueError` when any independent parameter is missing from the stored
data, and also verify that any stored *dependent* parameters match the values
recalculated from the independent parameters.

Log parsing
-----------

:func:`parse_micro_log` and :func:`parse_macro_log` parse Fortran simulator
output files (``micro_*.txt``, ``macro_*.txt``) and return a dict of
``{python_name: value}``.  Any numeric-valued Fortran name that cannot be
mapped to a Python parameter raises :exc:`ValueError`.

Log verification
----------------

:func:`verify_micro_params` and :func:`verify_macro_params` combine parsing
with comparison against a parameter object, raising :exc:`ValueError` on any
mismatch.

Aliases and overrides
---------------------

The public functions accept two optional keyword arguments:

* ``aliases`` — a ``{python_name: fortran_name}`` dict that tells the parser
  to look up an alternative Fortran name in the log file.  This also suppresses
  the "unknown Fortran name" error for the aliased entry.  Used by the parse
  and verify functions.
* ``overrides`` — a ``{python_name: value}`` dict of direct values (int, float,
  Quantity-string, etc.) used as-is.  Used by the load and verify functions.

Note on ``kon`` in macro logs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Fortran macroscale code writes ``kon=`` for the tPA binding rate, while the
microscale code writes ``ktPAon=``.  Only ``ktPAon`` appears in
:meth:`~.parameters.MicroParameters.fortran_names`.  When parsing macro logs
that contain ``kon=``, pass
``aliases={"bind_rate_tPA": "kon"}`` to resolve it.
"""

import dataclasses
import inspect
import re
from pathlib import Path

from pint import Quantity

from .constants import Q_
from .parameters import MacroParameters, MicroParameters


__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _class_fortran_names(cls):
    """Return ``{python_name: fortran_spec}`` restricted to fields of *cls*.

    :meth:`~.parameters.Parameters.fortran_names` parses the entire source
    file and therefore returns entries for both ``MicroParameters`` and
    ``MacroParameters``.  This helper filters to only the fields declared in
    *cls*.

    :param cls: A :class:`~.parameters.Parameters` subclass.
    :return: Filtered ``{python_name: fortran_spec}`` dict.
    :rtype: dict[str, str]
    """
    all_names = cls.fortran_names()
    cls_field_names = {f.name for f in dataclasses.fields(cls)}
    return {k: v for k, v in all_names.items() if k in cls_field_names}


def _build_inverse_map(cls, extra_cls=None):
    """Build ``{fortran_name_lower: (python_name, transform, source_cls)}``.

    Iterates :func:`_class_fortran_names` for *cls* (and optionally
    *extra_cls*) and inverts the mapping.  The primary class (*cls*) takes
    precedence: if both classes share the same Fortran name (e.g. ``seed``
    for both ``micro_seed`` and ``macro_seed``), the *cls* entry wins and the
    *extra_cls* entry is silently dropped.

    Transform labels encode the Fortran → Python conversion direction:

    * ``None``      — identity: ``python = fortran``
    * ``'minus1'``  — subtract: ``python = fortran − 1``
    * ``'times100'``— divide:   ``python = fortran / 100``

    :param cls: Primary parameter class.
    :param extra_cls: Optional additional class whose non-conflicting names
        are also included (for cross-class params in macro logs).
    :return: Inverse Fortran-name map.
    :rtype: dict[str, tuple]
    """
    inverse = {}

    def _add(c, overwrite):
        for py_name, fortran_spec in _class_fortran_names(c).items():
            if fortran_spec.endswith('-1'):
                base_name = fortran_spec[:-2]
                transform = 'minus1'
            elif fortran_spec.endswith('*100'):
                base_name = fortran_spec[:-4]
                transform = 'times100'
            else:
                base_name = fortran_spec
                transform = None
            key = base_name.lower()
            if overwrite or key not in inverse:
                inverse[key] = (py_name, transform, c)

    _add(cls, overwrite=True)
    if extra_cls is not None:
        _add(extra_cls, overwrite=False)

    return inverse


def _parse_fortran_kv(lines):
    """Extract ``key = value`` pairs from Fortran log lines.

    Matches two formats:

    * ``key = value`` — standard Fortran output (pattern
      ``r'^\\s*(\\w+)\\s*=\\s*(.+?)\\s*$'``).
    * ``Setting key = value`` — command-line argument echoes written by the
      Fortran code during startup (pattern
      ``r'^\\s*Setting\\s+(\\w+)\\s*=\\s*(.+?)\\s*$'``).

    Lines that match neither pattern (e.g. file-path lines, multi-line values)
    are silently skipped.  When the same key appears in both a ``Setting`` line
    and a later ``key = value`` line, the later value overwrites the earlier
    one.

    :param lines: Iterable of text lines from a Fortran log file.
    :return: ``{fortran_name_lower: raw_value_str}``.
    :rtype: dict[str, str]
    """
    kv_pattern = re.compile(r'^\s*(\w+)\s*=\s*(.+?)\s*$')
    setting_pattern = re.compile(r'^\s*Setting\s+(\w+)\s*=\s*(.+?)\s*$')
    result = {}
    for line in lines:
        m = setting_pattern.match(line) or kv_pattern.match(line)
        if m:
            result[m.group(1).lower()] = m.group(2).strip()
    return result



def _apply_transform(raw_num, transform):
    """Apply the Fortran → Python transform to a numeric value.

    :param raw_num: Float value read directly from the Fortran log.
    :param transform: ``None``, ``'minus1'``, or ``'times100'``.
    :return: Transformed numeric value.
    :rtype: float
    """
    if transform == 'minus1':
        return raw_num - 1
    if transform == 'times100':
        return raw_num / 100
    return raw_num


def _to_quantity_or_number(val, py_name, units_dict):
    """Wrap *val* in a :class:`~pint.Quantity` if *py_name* has units.

    :param val: Numeric value (post-transform).
    :param py_name: Python parameter name.
    :param units_dict: ``{py_name: unit_str}`` from ``cls.units()``.
    :return: ``Quantity`` or bare ``float``.
    """
    unit = units_dict.get(py_name)
    if unit is not None:
        return Q_(val, unit)
    return val


def _extract_numeric(value, units_dict, py_name):
    """Extract a comparable float magnitude from *value*.

    Handles :class:`~pint.Quantity` (converted to canonical units), plain
    ``int``/``float``, and strings (parsed with Pint if they contain units,
    otherwise ``float()``).  Returns ``None`` for types that cannot be
    meaningfully compared (``list``, ``tuple``, ``dict``, non-numeric
    strings).

    :param value: A stored or calculated parameter value.
    :param units_dict: ``{py_name: unit_str}`` for unit conversion.
    :param py_name: Python parameter name.
    :return: Float magnitude in canonical units, or ``None``.
    :rtype: float | None
    """
    unit = units_dict.get(py_name)

    if isinstance(value, Quantity):
        try:
            if unit:
                return float(value.to(unit).magnitude)
            return float(value.magnitude)
        except Exception:
            return None

    if isinstance(value, str):
        try:
            q = Q_(value)
            if unit:
                return float(q.to(unit).magnitude)
            return float(q.magnitude)
        except Exception:
            try:
                return float(value)
            except (ValueError, TypeError):
                return None

    if isinstance(value, (int, float)):
        return float(value)

    # list, tuple, dict, etc. — skip
    return None


def _compare_numeric(stored_num, calc_num, tolerance):
    """Return ``True`` if *stored_num* ≈ *calc_num* within *tolerance*.

    Uses relative tolerance unless *calc_num* is zero, in which case absolute
    tolerance is applied.

    :param stored_num: The stored (or log-parsed) numeric value.
    :param calc_num: The expected (calculated) numeric value.
    :param tolerance: Tolerance threshold.
    :rtype: bool
    """
    if calc_num == 0:
        return abs(stored_num - calc_num) < tolerance
    return abs(stored_num - calc_num) / abs(calc_num) < tolerance


def _check_dependent_params(base_params, instance, cls, tolerance):
    """Compare stored dependent parameters against recalculated values.

    Iterates over keys in *base_params* that are **not** in the independent
    parameter set for *cls*.  Skips complex types (``list``, ``tuple``,
    ``dict``, ``str``) and any value that cannot be extracted as a number.

    :param base_params: Raw dict as loaded from HDF5 (base types only).
    :param instance: Freshly built parameter dataclass.
    :param cls: The parameter class (``MicroParameters`` or
        ``MacroParameters``).
    :param tolerance: Numeric tolerance.
    :return: List of ``(name, stored_value, calculated_value)`` for each
        mismatch found.
    :rtype: list[tuple]
    """
    independent = set(inspect.signature(cls).parameters.keys())
    units_dict = cls.units()
    mismatches = []

    for name, stored_val in base_params.items():
        if name in independent:
            continue
        if not hasattr(instance, name):
            continue
        calculated_val = getattr(instance, name)
        # Skip complex types that cannot be numerically compared
        if isinstance(calculated_val, (list, tuple, dict, str)):
            continue

        stored_num = _extract_numeric(stored_val, units_dict, name)
        calc_num = _extract_numeric(calculated_val, units_dict, name)

        if stored_num is None or calc_num is None:
            continue
        if not _compare_numeric(stored_num, calc_num, tolerance):
            mismatches.append((name, stored_val, calculated_val))

    return mismatches


# ─────────────────────────────────────────────────────────────────────────────
# Strict loading functions
# ─────────────────────────────────────────────────────────────────────────────

def load_micro_params(base_params, *, overrides=None, tolerance=1e-9):
    """Strictly load and validate :class:`~.parameters.MicroParameters`.

    Unlike :meth:`~.parameters.Parameters.parse_from_basedict`, this function:

    * Raises :exc:`ValueError` if any independent parameter is missing.
    * Verifies that any stored dependent parameters match those recalculated
      from the independent parameters.

    :param base_params: Dict of parameter values (base Python types only), as
        returned by :meth:`~.parameters.Parameters.to_basedict` or loaded from
        an HDF5 attribute group.
    :type base_params: dict[str, int | float | str]
    :param overrides: Optional ``{python_name: value}`` to supply missing
        independent parameters or replace stored values.  Values are used
        directly (int, float, Quantity-string, etc.).
    :type overrides: dict | None
    :param tolerance: Relative tolerance for comparing dependent parameter
        values (absolute when the calculated value is zero).
    :type tolerance: float
    :return: Validated :class:`~.parameters.MicroParameters` instance.
    :rtype: MicroParameters
    :raises ValueError: If required independent parameters are missing, or if
        stored dependent parameters do not match recalculated values.
    """
    direct_overrides = overrides if overrides else {}

    independent = set(inspect.signature(MicroParameters).parameters.keys())

    # Start from stored values for independent params, then apply direct overrides
    merged = {k: v for k, v in base_params.items() if k in independent}
    for k, v in direct_overrides.items():
        if k in independent:
            merged[k] = v

    missing = independent - set(merged.keys())
    if missing:
        raise ValueError(
            f"Missing independent MicroParameters: {sorted(missing)}"
        )

    instance = MicroParameters.parse_from_basedict(merged)

    mismatches = _check_dependent_params(base_params, instance, MicroParameters, tolerance)
    if mismatches:
        details = "; ".join(
            f"{name}: stored {stored!r} != calculated {calc!r}"
            for name, stored, calc in mismatches
        )
        raise ValueError(
            f"Dependent parameter mismatch in MicroParameters: {details}"
        )

    return instance


def load_macro_params(base_params, micro_params, *, overrides=None, tolerance=1e-9):
    """Strictly load and validate :class:`~.parameters.MacroParameters`.

    Like :func:`load_micro_params` but for macroscale parameters.  The
    *micro_params* argument provides the associated microscale parameter set,
    which is required to construct a :class:`~.parameters.MacroParameters`
    instance but is not stored in *base_params*.

    :param base_params: Dict of macroscale parameter values (base Python types
        only).
    :type base_params: dict[str, int | float | str]
    :param micro_params: The associated :class:`~.parameters.MicroParameters`
        instance.  Must not be ``None``.
    :type micro_params: MicroParameters
    :param overrides: Optional ``{python_name: value}`` overrides.  Values are
        used directly (int, float, Quantity-string, etc.).
    :type overrides: dict | None
    :param tolerance: Relative tolerance for dependent parameter comparison.
    :type tolerance: float
    :return: Validated :class:`~.parameters.MacroParameters` instance.
    :rtype: MacroParameters
    :raises ValueError: If *micro_params* is ``None``, if required independent
        parameters are missing, or if stored dependent parameters do not match.
    """
    if micro_params is None:
        raise ValueError("micro_params is required for load_macro_params")

    direct_overrides = overrides if overrides else {}

    independent = set(inspect.signature(MacroParameters).parameters.keys())
    independent.discard('micro_params')  # Injected separately, not in base_params

    merged = {k: v for k, v in base_params.items() if k in independent}
    for k, v in direct_overrides.items():
        if k in independent:
            merged[k] = v

    missing = independent - set(merged.keys())
    if missing:
        raise ValueError(
            f"Missing independent MacroParameters: {sorted(missing)}"
        )

    # Inject micro_params so parse_from_basedict can construct the instance
    merged['micro_params'] = micro_params
    instance = MacroParameters.parse_from_basedict(merged)

    mismatches = _check_dependent_params(base_params, instance, MacroParameters, tolerance)
    if mismatches:
        details = "; ".join(
            f"{name}: stored {stored!r} != calculated {calc!r}"
            for name, stored, calc in mismatches
        )
        raise ValueError(
            f"Dependent parameter mismatch in MacroParameters: {details}"
        )

    return instance


# ─────────────────────────────────────────────────────────────────────────────
# Internal log-parsing core
# ─────────────────────────────────────────────────────────────────────────────

def _parse_log(path, stop_pattern, inverse_map, units_dict, alias_map):
    """Core Fortran log parser shared by micro and macro parsers.

    Reads *path* line by line, stopping before the first line matching
    *stop_pattern*.  Extracts ``key = value`` pairs, maps them via
    *inverse_map* (augmented with *alias_map*), and applies the appropriate
    Fortran → Python transform.  Raises :exc:`ValueError` for any
    numeric-valued Fortran name that cannot be resolved.

    Non-numeric values (e.g. ``filetype=binary``) are silently skipped
    regardless of whether the key is in the inverse map.

    :param path: Path to the Fortran log file.
    :type path: str | Path
    :param stop_pattern: Compiled :mod:`re` pattern; parsing stops at the
        first line that matches.
    :param inverse_map: ``{fortran_lower: (py_name, transform, source_cls)}``.
    :param units_dict: Combined ``{py_name: unit_str}`` for unit wrapping.
    :param alias_map: ``{fortran_alias_lower: py_name}`` from caller overrides.
    :return: ``{python_name: value}`` for all resolved parameters.
    :rtype: dict
    :raises ValueError: If any line contains a numeric-valued Fortran name
        that is not in *inverse_map* and not covered by *alias_map*.
    """
    path = Path(path)
    with path.open('r', errors='replace') as fh:
        lines = []
        for line in fh:
            if stop_pattern.search(line):
                break
            lines.append(line)

    raw = _parse_fortran_kv(lines)

    result = {}
    unknown = []

    for fortran_lower, raw_str in raw.items():
        # Alias overrides take priority
        if fortran_lower in alias_map:
            py_name = alias_map[fortran_lower]
            try:
                num = float(raw_str)
            except (ValueError, TypeError):
                continue  # Non-numeric; skip silently
            result[py_name] = _to_quantity_or_number(num, py_name, units_dict)
            continue

        if fortran_lower in inverse_map:
            py_name, transform, _ = inverse_map[fortran_lower]
            try:
                num = float(raw_str)
            except (ValueError, TypeError):
                continue  # Non-numeric; skip silently (e.g. filetype=binary)
            transformed = _apply_transform(num, transform)
            result[py_name] = _to_quantity_or_number(transformed, py_name, units_dict)
            continue

        # Unknown key — flag only if the value is numeric (a real parameter)
        try:
            float(raw_str)
            unknown.append(fortran_lower)
        except (ValueError, TypeError):
            pass  # Non-numeric unknown keys are silently ignored

    if unknown:
        raise ValueError(
            f"Unrecognized Fortran parameter names in log: {sorted(unknown)}"
        )

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Public log-parsing functions
# ─────────────────────────────────────────────────────────────────────────────

def parse_micro_log(path, *, aliases=None):
    """Parse a Fortran microscale log file into Python parameter values.

    Reads from the start of *path* until the first ``stats =`` line (exclusive).
    Maps Fortran variable names to Python parameter names using
    :meth:`~.parameters.MicroParameters.fortran_names` metadata.

    :param path: Path to the micro log file
        (e.g. ``micro_PLG2_tPA01_TB-xiii.txt``).
    :type path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases.
        Each entry causes *fortran_name* in the log to map to *python_name*
        instead of raising an unknown-name error.
    :type aliases: dict | None
    :return: ``{python_name: value}`` for all recognized parameters.
    :rtype: dict
    :raises ValueError: If any numeric-valued Fortran name cannot be mapped to
        a Python parameter.
    """
    alias_map = {v.lower(): k for k, v in aliases.items()} if aliases else {}
    inverse_map = _build_inverse_map(MicroParameters)
    units_dict = MicroParameters.units()
    stop_pattern = re.compile(r'^\s*stats\s*=')
    return _parse_log(path, stop_pattern, inverse_map, units_dict, alias_map)


def parse_macro_log(path, *, aliases=None):
    """Parse a Fortran macroscale log file into Python parameter values.

    Reads from the start of *path* until the first line beginning with
    ``After `` (exclusive).  Maps Fortran variable names using both
    :meth:`~.parameters.MacroParameters.fortran_names` and
    :meth:`~.parameters.MicroParameters.fortran_names`.
    Cross-class parameters that appear in the macro log (e.g. ``bs``) but
    belong to :class:`~.parameters.MicroParameters` are included in the result
    dict; callers can identify them by checking against
    ``dataclasses.fields(MicroParameters)``.

    .. note::
        The macro Fortran code writes ``kon=`` for the tPA binding rate, while
        the microscale code uses ``ktPAon=``.  The parameter
        ``bind_rate_tPA`` maps to ``ktPAon`` in
        :meth:`~.parameters.MicroParameters.fortran_names`, so ``kon=`` lines
        in macro logs will raise :exc:`ValueError` unless you pass
        ``aliases={"bind_rate_tPA": "kon"}``.

    :param path: Path to the macro log file
        (e.g. ``macro_TB-xiii__21_105_00.txt``).
    :type path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases.
    :type aliases: dict | None
    :return: ``{python_name: value}`` for all recognized parameters.
    :rtype: dict
    :raises ValueError: If any numeric-valued Fortran name cannot be mapped.
    """
    alias_map = {v.lower(): k for k, v in aliases.items()} if aliases else {}
    # MacroParameters takes precedence for shared names (e.g. seed, simulations)
    inverse_map = _build_inverse_map(MacroParameters, extra_cls=MicroParameters)
    units_dict = {**MicroParameters.units(), **MacroParameters.units()}
    stop_pattern = re.compile(r'^After\s')
    return _parse_log(path, stop_pattern, inverse_map, units_dict, alias_map)


# ─────────────────────────────────────────────────────────────────────────────
# Internal verification core
# ─────────────────────────────────────────────────────────────────────────────

def _verify_params(params_obj, log_params, units_dict, tolerance, direct_overrides):
    """Core verification logic comparing log values against a parameter object.

    :param params_obj: The parameter dataclass instance to check against.
    :param log_params: ``{python_name: value}`` parsed from the Fortran log.
    :param units_dict: Combined ``{py_name: unit_str}`` for numeric extraction.
    :param tolerance: Numeric tolerance.
    :param direct_overrides: ``{py_name: value}`` direct overrides from the
        caller; these replace the stored value for comparison purposes.
    :return: List of ``(name, log_value, stored_value)`` mismatches.
    :rtype: list[tuple]
    """
    mismatches = []

    for py_name, log_val in log_params.items():
        if py_name in direct_overrides:
            comparison_val = direct_overrides[py_name]
        elif hasattr(params_obj, py_name):
            comparison_val = getattr(params_obj, py_name)
        else:
            continue

        log_num = _extract_numeric(log_val, units_dict, py_name)
        stored_num = _extract_numeric(comparison_val, units_dict, py_name)

        if log_num is None or stored_num is None:
            continue
        if not _compare_numeric(log_num, stored_num, tolerance):
            mismatches.append((py_name, log_val, comparison_val))

    return mismatches


# ─────────────────────────────────────────────────────────────────────────────
# Public verification functions
# ─────────────────────────────────────────────────────────────────────────────

def verify_micro_params(micro_params, log_path, *, aliases=None, overrides=None, tolerance=1e-9):
    """Verify that a :class:`~.parameters.MicroParameters` matches a micro log.

    Parses *log_path* and compares every recognized parameter value against
    *micro_params*.  All mismatches are collected and reported together in a
    single :exc:`ValueError`.

    :param micro_params: The parameter instance to validate.
    :type micro_params: MicroParameters
    :param log_path: Path to the Fortran micro log file.
    :type log_path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases for the
        log parser.  Each entry causes *fortran_name* in the log to map to
        *python_name* instead of raising an unknown-name error.
    :type aliases: dict | None
    :param overrides: Optional ``{python_name: value}`` direct overrides.
        Each entry replaces the stored attribute in the comparison (e.g. to
        accept a known corrected value).
    :type overrides: dict | None
    :param tolerance: Relative (or absolute when zero) numeric tolerance.
    :type tolerance: float
    :return: ``None`` on success.
    :raises ValueError: If any log parameter value does not match
        *micro_params* within *tolerance*.
    """
    direct_overrides = overrides if overrides else {}
    log_params = parse_micro_log(log_path, aliases=aliases)
    units_dict = MicroParameters.units()

    mismatches = _verify_params(
        micro_params, log_params, units_dict, tolerance, direct_overrides
    )
    if mismatches:
        details = "; ".join(
            f"{name}: log {log_val!r} != stored {stored!r}"
            for name, log_val, stored in mismatches
        )
        raise ValueError(f"MicroParameters do not match log: {details}")


def verify_macro_params(macro_params, log_path, *, aliases=None, overrides=None, tolerance=1e-9):
    """Verify that a :class:`~.parameters.MacroParameters` matches a macro log.

    Like :func:`verify_micro_params` but for macroscale parameters.
    Cross-class parameters that belong to
    :class:`~.parameters.MicroParameters` (e.g. ``bs``) are verified against
    ``macro_params.micro_params``.

    :param macro_params: The parameter instance to validate.
    :type macro_params: MacroParameters
    :param log_path: Path to the Fortran macro log file.
    :type log_path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases for the
        log parser.
    :type aliases: dict | None
    :param overrides: Optional ``{python_name: value}`` direct overrides.
    :type overrides: dict | None
    :param tolerance: Relative (or absolute when zero) numeric tolerance.
    :type tolerance: float
    :return: ``None`` on success.
    :raises ValueError: If any log parameter value does not match
        *macro_params* (or its associated *micro_params*) within *tolerance*.
    """
    direct_overrides = overrides if overrides else {}
    log_params = parse_macro_log(log_path, aliases=aliases)

    units_dict = {**MicroParameters.units(), **MacroParameters.units()}
    micro_field_names = {f.name for f in dataclasses.fields(MicroParameters)}

    # Partition log_params: cross-class (MicroParameters) vs MacroParameters
    macro_log = {}
    micro_log = {}
    for py_name, val in log_params.items():
        if py_name in micro_field_names:
            micro_log[py_name] = val
        else:
            macro_log[py_name] = val

    mismatches = _verify_params(
        macro_params, macro_log, units_dict, tolerance, direct_overrides
    )

    # Cross-class params are verified against macro_params.micro_params
    if micro_log and macro_params.micro_params is not None:
        micro_mismatches = _verify_params(
            macro_params.micro_params,
            micro_log,
            MicroParameters.units(),
            tolerance,
            direct_overrides,
        )
        mismatches.extend(micro_mismatches)

    if mismatches:
        details = "; ".join(
            f"{name}: log {log_val!r} != stored {stored!r}"
            for name, log_val, stored in mismatches
        )
        raise ValueError(f"MacroParameters do not match log: {details}")
