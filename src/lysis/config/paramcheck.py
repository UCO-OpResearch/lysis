"""Parameter validation and Fortran log verification.

Provides strict parameter loading and Fortran log verification on top
of :meth:`~.parameters.Parameters.parse_from_basedict`.  All functions here
are **optional** — existing code that calls ``parse_from_basedict()`` directly
is unchanged.

Strict loading
--------------

:func:`load_micro_params` and :func:`load_macro_params` raise
:exc:`ValueError` when any independent parameter is missing from the stored
data, and also verify that any stored *dependent* parameters match the values
recalculated from the independent parameters.

Log verification
----------------

:func:`verify_micro_params` and :func:`verify_macro_params` parse Fortran
log files (via :mod:`lysis.dataio.fileops`), resolve Fortran names to Python
names using :meth:`~.parameters.Parameters.inverse_fortran_map`, and compare
the values against a parameter object.  All mismatches are collected and
reported together in a single :exc:`ValueError`.

Aliases and overrides
---------------------

* ``aliases`` — a ``{python_name: fortran_name}`` dict that tells the
  verify functions to look up an alternative Fortran name in the parsed
  log output.  Used by :func:`verify_micro_params` and
  :func:`verify_macro_params`.
* ``overrides`` — a ``{python_name: value}`` dict of direct values (int,
  float, Quantity-string, etc.) used as-is.  Used by the load and verify
  functions.
"""

import dataclasses
import inspect

from pint import Quantity

from .constants import Q_
from .parameters import MacroParameters, MicroParameters, Parameters


__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


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

    * Raises :exc:`ValueError` if any independent parameter is missing
      *and has no default value*.
    * Verifies that any stored dependent parameters match those recalculated
      from the independent parameters.

    Parameters that have default values in :class:`~.parameters.MicroParameters`
    but are not present in *base_params* or *overrides* are filled in
    automatically.  This covers parameters with no Fortran equivalent
    (e.g. ``fibrinogen_length``, ``micro_version``) that never appear in
    Fortran log files.

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

    # Fill in defaults for missing parameters that have default values.
    # This covers parameters with no Fortran equivalent (e.g.
    # fibrinogen_length, micro_version) that never appear in log files.
    sig = inspect.signature(MicroParameters)
    units_dict = MicroParameters.units()
    for param_name, param in sig.parameters.items():
        if param_name not in merged and param.default is not inspect.Parameter.empty:
            default = param.default
            if isinstance(default, Quantity):
                unit = units_dict.get(param_name)
                if unit:
                    merged[param_name] = str(default.to(unit))
                else:
                    merged[param_name] = str(default)
            else:
                merged[param_name] = default

    missing = independent - set(merged.keys())
    if missing:
        raise ValueError(f"Missing independent MicroParameters: {sorted(missing)}")

    instance = MicroParameters.parse_from_basedict(merged)

    mismatches = _check_dependent_params(
        base_params, instance, MicroParameters, tolerance
    )
    if mismatches:
        details = "; ".join(
            f"{name}: stored {stored!r} != calculated {calc!r}"
            for name, stored, calc in mismatches
        )
        raise ValueError(f"Dependent parameter mismatch in MicroParameters: {details}")

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
    independent.discard("micro_params")  # Injected separately, not in base_params

    merged = {k: v for k, v in base_params.items() if k in independent}
    for k, v in direct_overrides.items():
        if k in independent:
            merged[k] = v

    # Fill in defaults for missing parameters (same rationale as load_micro_params)
    sig = inspect.signature(MacroParameters)
    units_dict = MacroParameters.units()
    for param_name, param in sig.parameters.items():
        if param_name == "micro_params":
            continue
        if param_name not in merged and param.default is not inspect.Parameter.empty:
            default = param.default
            if isinstance(default, Quantity):
                unit = units_dict.get(param_name)
                if unit:
                    merged[param_name] = str(default.to(unit))
                else:
                    merged[param_name] = str(default)
            else:
                merged[param_name] = default

    missing = independent - set(merged.keys())
    if missing:
        raise ValueError(f"Missing independent MacroParameters: {sorted(missing)}")

    # Inject micro_params so parse_from_basedict can construct the instance
    merged["micro_params"] = micro_params
    instance = MacroParameters.parse_from_basedict(merged)

    mismatches = _check_dependent_params(
        base_params, instance, MacroParameters, tolerance
    )
    if mismatches:
        details = "; ".join(
            f"{name}: stored {stored!r} != calculated {calc!r}"
            for name, stored, calc in mismatches
        )
        raise ValueError(f"Dependent parameter mismatch in MacroParameters: {details}")

    return instance


# ─────────────────────────────────────────────────────────────────────────────
# Internal verification helpers
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_log_params(raw_log, inverse_map, units_dict, alias_map):
    """Resolve raw Fortran log params to ``{python_name: value}`` dict.

    :param raw_log: ``{fortran_name_lower: float}`` from
        :func:`~lysis.dataio.fileops.parse_micro_log` or
        :func:`~lysis.dataio.fileops.parse_macro_log`.
    :param inverse_map: From :meth:`Parameters.inverse_fortran_map`.
    :param units_dict: ``{py_name: unit_str}`` for unit wrapping.
    :param alias_map: ``{fortran_alias_lower: py_name}`` from caller.
    :return: ``{python_name: value}`` for all resolved parameters.
    :rtype: dict
    """
    result = {}

    for fortran_lower, raw_num in raw_log.items():
        # Alias overrides take priority
        if fortran_lower in alias_map:
            py_name = alias_map[fortran_lower]
            result[py_name] = Parameters.to_quantity_or_number(
                raw_num, py_name, units_dict
            )
            continue

        if fortran_lower in inverse_map:
            py_name, transform, _ = inverse_map[fortran_lower]
            transformed = Parameters.apply_fortran_transform(raw_num, transform)
            result[py_name] = Parameters.to_quantity_or_number(
                transformed, py_name, units_dict
            )
            continue

        # Unknown key — silently skip (raw log already filtered to numeric)

    return result


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


def verify_micro_params(
    micro_params, log_path, *, aliases=None, overrides=None, tolerance=1e-9
):
    """Verify that a :class:`~.parameters.MicroParameters` matches a micro log.

    Parses *log_path* via :func:`~lysis.dataio.fileops.parse_micro_log`,
    resolves Fortran names using :meth:`MicroParameters.inverse_fortran_map`,
    and compares every recognized parameter value against *micro_params*.
    All mismatches are collected and reported together in a single
    :exc:`ValueError`.

    :param micro_params: The parameter instance to validate.
    :type micro_params: MicroParameters
    :param log_path: Path to the Fortran micro log file.
    :type log_path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases.
        Each entry causes *fortran_name* in the log to map to *python_name*
        instead of being ignored.
    :type aliases: dict | None
    :param overrides: Optional ``{python_name: value}`` direct overrides.
        Each entry replaces the stored attribute in the comparison.
    :type overrides: dict | None
    :param tolerance: Relative (or absolute when zero) numeric tolerance.
    :type tolerance: float
    :return: ``None`` on success.
    :raises ValueError: If any log parameter value does not match
        *micro_params* within *tolerance*.
    """
    from ..dataio.fileops import parse_micro_log

    direct_overrides = overrides if overrides else {}
    alias_map = {v.lower(): k for k, v in aliases.items()} if aliases else {}

    raw_log = parse_micro_log(log_path)
    inverse_map = MicroParameters.inverse_fortran_map()
    units_dict = MicroParameters.units()

    log_params = _resolve_log_params(raw_log, inverse_map, units_dict, alias_map)

    mismatches = _verify_params(
        micro_params, log_params, units_dict, tolerance, direct_overrides
    )
    if mismatches:
        details = "; ".join(
            f"{name}: log {log_val!r} != stored {stored!r}"
            for name, log_val, stored in mismatches
        )
        raise ValueError(f"MicroParameters do not match log: {details}")


def verify_macro_params(
    macro_params, log_path, *, aliases=None, overrides=None, tolerance=1e-9
):
    """Verify that a :class:`~.parameters.MacroParameters` matches a macro log.

    Like :func:`verify_micro_params` but for macroscale parameters.
    Cross-class parameters that belong to
    :class:`~.parameters.MicroParameters` (e.g. ``bs``) are verified against
    ``macro_params.micro_params``.

    :param macro_params: The parameter instance to validate.
    :type macro_params: MacroParameters
    :param log_path: Path to the Fortran macro log file.
    :type log_path: str | Path
    :param aliases: Optional ``{python_name: fortran_name}`` aliases.
    :type aliases: dict | None
    :param overrides: Optional ``{python_name: value}`` direct overrides.
    :type overrides: dict | None
    :param tolerance: Relative (or absolute when zero) numeric tolerance.
    :type tolerance: float
    :return: ``None`` on success.
    :raises ValueError: If any log parameter value does not match
        *macro_params* (or its associated *micro_params*) within *tolerance*.
    """
    from ..dataio.fileops import parse_macro_log

    direct_overrides = overrides if overrides else {}
    alias_map = {v.lower(): k for k, v in aliases.items()} if aliases else {}

    raw_log = parse_macro_log(log_path)
    inverse_map = MacroParameters.inverse_fortran_map(extra_cls=MicroParameters)
    units_dict = {**MicroParameters.units(), **MacroParameters.units()}

    log_params = _resolve_log_params(raw_log, inverse_map, units_dict, alias_map)

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
