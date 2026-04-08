"""Algebraic parameter resolution for Experiment initialization.

Given a partial set of parameter values — which may include dependent
(computed) parameters in place of independent (input) ones — this module
solves for any algebraically determinable missing values.

Typical usage
-------------

::

    from lysis.config.param_resolver import resolve_micro_params, resolve_macro_params
    from lysis.config.paramcheck import load_micro_params, load_macro_params

    # Resolve algebraic relationships (e.g. solve for bind_rate from koff/Kd)
    resolved_micro = resolve_micro_params(csv_row_micro_dict)
    micro_params = load_micro_params(resolved_micro)

    resolved_macro = resolve_macro_params(csv_row_macro_dict, micro_params.to_basedict())
    macro_params = load_macro_params(resolved_macro, micro_params)

Design notes
------------

All algebraic equations are encoded with the same canonical units used by
:meth:`~.parameters.Parameters.to_basedict` (i.e. those listed in each
parameter's ``:Units:`` docstring tag).  Before being handed to ``sympy``,
every provided :class:`~pint.Quantity` is converted to its canonical unit
magnitude.  Solved magnitudes are returned as ``"{magnitude} {unit}"``
strings so they are directly consumable by
:meth:`~.parameters.Parameters.parse_from_basedict`.
"""

import inspect
import warnings
from typing import Any

import sympy
from pint import Quantity

from .constants import Q_, ureg
from .parameters import MacroParameters, MicroParameters, Parameters


__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


__all__ = [
    "resolve_micro_params",
    "resolve_macro_params",
    "ParameterConflict",
    "UnderdeterminedParameters",
]


# ─── Exceptions ──────────────────────────────────────────────────────────────


class ParameterConflict(ValueError):
    """Raised when provided parameter values are mutually inconsistent.

    :param conflicts: List of ``(param_name, provided_value, computed_value)``
        triples, one per conflicting parameter.
    :type conflicts: list[tuple[str, Any, Any]]
    """

    def __init__(self, conflicts: list[tuple[str, Any, Any]]):
        self.conflicts = conflicts
        lines = [
            f"  {name}: provided {provided!r}, computed {computed!r}"
            for name, provided, computed in conflicts
        ]
        super().__init__("Parameter value conflicts:\n" + "\n".join(lines))


class UnderdeterminedParameters(ValueError):
    """Raised when required parameters cannot be determined from the provided values.

    :param unresolved: List of parameter names that could not be resolved.
    :type unresolved: list[str]
    """

    def __init__(self, unresolved: list[str]):
        self.unresolved = sorted(unresolved)
        super().__init__(
            f"Cannot resolve parameters (insufficient constraints): {self.unresolved}"
        )


# ─── Physical constants ───────────────────────────────────────────────────────

# Avogadro's constant in mol^-1, taken from Pint's unit registry to ensure
# exact agreement with MicroParameters.__post_init__ calculations.
_N_A = float(Q_(1, ureg.avogadro_constant).to("1/mol").magnitude)

# Conversion factor from mol/µm³ to µM.
# Derivation: 1 µm = 1e-6 m; 1 µm³ = 1e-18 m³ = 1e-15 L
# → 1 mol/µm³ = 1/(1e-15 L) mol = 1e15 M = 1e21 µM
_MOL_PER_UM3_TO_MICROMOLAR = 1e21


# ─── Canonical units ─────────────────────────────────────────────────────────

# Cache the units dict once at module load (regex parse of parameters.py).
_UNITS: dict[str, str] = Parameters.units()


# ─── Symbol tables ────────────────────────────────────────────────────────────

# Parameters that appear in the MicroParameters algebraic equations.
# This includes both independent (init=True) and dependent (init=False) params
# that have algebraic relationships worth encoding.
_MICRO_EQ_PARAMS: list[str] = [
    "fibrinogen_length",
    "fibrinogen_radius",
    "fiber_radius",
    "protofibril_radius",
    "nodes_in_micro_row",
    "diss_const_tPA_wPLG",
    "diss_const_tPA_woPLG",
    "diss_const_PLG_intact",
    "diss_const_PLG_nicked",
    "bind_rate_tPA",
    "bind_rate_PLG",
    "unbind_rate_PLG_intact",
    "unbind_rate_PLG_nicked",
    "unbind_rate_tPA_wPLG",
    "unbind_rate_tPA_woPLG",
    "protein_per_fiber",
    "fibrin_conc_per_fiber",
    "binding_sites",
]

# Parameters that appear in the MacroParameters algebraic equations.
_MACRO_EQ_PARAMS: list[str] = [
    "cols",
    "rows",
    "empty_rows",
    "full_row",
    "xz_row",
    "fiber_rows",
    "empty_edges",
    "total_edges",
    "total_fibers",
    "moving_probability",
    "pore_size",
    "diffusion_coeff",
    "time_step",
]

_MICRO_SYMS: dict[str, sympy.Symbol] = {
    name: sympy.Symbol(name) for name in _MICRO_EQ_PARAMS
}
_MACRO_SYMS: dict[str, sympy.Symbol] = {
    name: sympy.Symbol(name) for name in _MACRO_EQ_PARAMS
}


# ─── Equation factories ───────────────────────────────────────────────────────


def _micro_equations() -> list[sympy.Eq]:
    """Return the 8 algebraic equations for MicroParameters dependent params.

    All values are expressed in canonical units (from ``Parameters.units()``):

    - Lengths: microns
    - Concentrations: micromolar
    - Binding/unbinding rates: sec⁻¹ or (micromolar·sec)⁻¹
    - ``protein_per_fiber``: %
    """
    s = _MICRO_SYMS
    return [
        # protofibril_radius [µm] = 2 × fibrinogen_radius [µm]
        sympy.Eq(s["protofibril_radius"], 2 * s["fibrinogen_radius"]),
        # Dissociation constant = unbinding rate / binding rate  ⟹  koff = kon × Kd
        sympy.Eq(
            s["unbind_rate_PLG_intact"],
            s["bind_rate_PLG"] * s["diss_const_PLG_intact"],
        ),
        sympy.Eq(
            s["unbind_rate_PLG_nicked"],
            s["bind_rate_PLG"] * s["diss_const_PLG_nicked"],
        ),
        sympy.Eq(
            s["unbind_rate_tPA_wPLG"],
            s["bind_rate_tPA"] * s["diss_const_tPA_wPLG"],
        ),
        sympy.Eq(
            s["unbind_rate_tPA_woPLG"],
            s["bind_rate_tPA"] * s["diss_const_tPA_woPLG"],
        ),
        # Volume fraction of protein per fiber [%].
        # Equation S1 from Bannish et al. 2017, simplified:
        #   nodes² × protofibril_radius² / fiber_radius² × 100
        # (fibrinogen_length/2 cancels; multiply by 100 to convert to percent)
        sympy.Eq(
            s["protein_per_fiber"],
            s["nodes_in_micro_row"] ** 2
            * s["protofibril_radius"] ** 2
            / s["fiber_radius"] ** 2
            * 100,
        ),
        # Fibrin concentration per fiber [µM].
        # Equation S2 from Bannish et al. 2017.
        # All lengths in µm; multiply by _MOL_PER_UM3_TO_MICROMOLAR to reach µM.
        sympy.Eq(
            s["fibrin_conc_per_fiber"],
            s["nodes_in_micro_row"] ** 2
            / (
                s["fibrinogen_length"] / 2
                * sympy.pi
                * s["fiber_radius"] ** 2
                * _N_A
            )
            * _MOL_PER_UM3_TO_MICROMOLAR,
        ),
        # Binding site concentration [µM].
        # Equation S3 from Bannish et al. 2017.
        sympy.Eq(
            s["binding_sites"],
            4
            * (s["nodes_in_micro_row"] - 1)
            / s["nodes_in_micro_row"] ** 2
            * s["fibrin_conc_per_fiber"],
        ),
    ]


def _macro_equations() -> list[sympy.Eq]:
    """Return the algebraic equations for MacroParameters dependent params.

    Grid geometry params are dimensionless integers.
    ``time_step`` uses pore_size [cm], diffusion_coeff [cm²/s] → time_step [s].
    """
    s = _MACRO_SYMS
    return [
        # Grid geometry
        sympy.Eq(s["full_row"], 3 * s["cols"] - 1),
        sympy.Eq(s["xz_row"], 2 * s["cols"] - 1),
        sympy.Eq(s["fiber_rows"], s["rows"] - s["empty_rows"]),
        sympy.Eq(s["empty_edges"], s["full_row"] * s["empty_rows"]),
        sympy.Eq(
            s["total_edges"],
            s["full_row"] * (s["rows"] - 1) + s["xz_row"],
        ),
        sympy.Eq(
            s["total_fibers"],
            s["full_row"] * (s["rows"] - s["empty_rows"] - 1) + s["xz_row"],
        ),
        # Timestep from Bannish et al. 2014, Equation 2.4
        # time_step [s] = moving_probability × pore_size[cm]² / (12 × diffusion_coeff[cm²/s])
        sympy.Eq(
            s["time_step"],
            s["moving_probability"] * s["pore_size"] ** 2 / (12 * s["diffusion_coeff"]),
        ),
    ]


# ─── Internal helpers ────────────────────────────────────────────────────────


def _parse_to_magnitude(name: str, value: Any) -> float:
    """Convert *value* to a float magnitude in the canonical unit for *name*.

    :param name: Python parameter name (used to look up canonical unit).
    :param value: The parameter value as a :class:`~pint.Quantity`, a
        unit-bearing string (e.g. ``"72.7 nm"``), or a plain number.
    :returns: Float magnitude in canonical units.
    :raises ValueError: If the value cannot be parsed or converted.
    """
    unit = _UNITS.get(name)

    if isinstance(value, Quantity):
        if unit:
            return float(value.to(unit).magnitude)
        return float(value.magnitude)

    if isinstance(value, str):
        try:
            q = Q_(value)
            if unit:
                return float(q.to(unit).magnitude)
            return float(q.magnitude)
        except Exception:
            # No unit in string; assume canonical unit with a warning
            try:
                mag = float(value)
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    f"Cannot parse '{name}' = {value!r} as a number or Quantity"
                ) from exc
            if unit:
                warnings.warn(
                    f"Parameter '{name}' has no units in provided value {value!r}; "
                    f"assuming canonical unit '{unit}'.",
                    RuntimeWarning,
                    stacklevel=4,
                )
            return mag

    if isinstance(value, (int, float)):
        if unit:
            warnings.warn(
                f"Parameter '{name}' = {value!r} has no units; "
                f"assuming canonical unit '{unit}'.",
                RuntimeWarning,
                stacklevel=4,
            )
        return float(value)

    raise ValueError(
        f"Cannot convert '{name}' = {value!r} (type {type(value).__name__}) to float"
    )


def _mag_to_str(name: str, magnitude: float) -> str | float | int:
    """Format *magnitude* as a Pint-compatible string in canonical units.

    Returns the magnitude unchanged (as ``int`` or ``float``) for dimensionless
    parameters.

    :param name: Python parameter name.
    :param magnitude: Numeric magnitude in canonical units.
    """
    unit = _UNITS.get(name)
    if unit is None:
        # Dimensionless: preserve integer type for int-typed params
        int_mag = int(round(magnitude))
        if abs(int_mag - magnitude) < 1e-9:
            return int_mag
        return magnitude
    return f"{magnitude} {unit}"


def _resolve_equations(
    equations: list[sympy.Eq],
    syms: dict[str, sympy.Symbol],
    eq_param_names: list[str],
    provided: dict[str, Any],
    tolerance: float = 1e-9,
) -> dict[str, float]:
    """Substitute provided values into *equations* and solve for unknowns.

    :param equations: Sympy equations relating the parameters.
    :param syms: ``{name: Symbol}`` dict covering all params in equations.
    :param eq_param_names: Names of all params that appear in equations.
    :param provided: Raw provided dict (any value types); only params in
        *eq_param_names* are used.
    :param tolerance: Relative tolerance for conflict detection.
    :returns: ``{name: magnitude}`` for all params that could be determined
        (provided + solved).
    :raises ParameterConflict: If fully-determined equations are inconsistent.
    """
    # 1. Parse provided values for equation params → canonical magnitudes
    provided_mags: dict[str, float] = {}
    for name in eq_param_names:
        val = provided.get(name)
        if val is None or val == "":
            continue
        try:
            provided_mags[name] = _parse_to_magnitude(name, val)
        except Exception as exc:
            raise ValueError(f"Error parsing parameter '{name}': {exc}") from exc

    # 2. Substitute all known values into equations
    subs = {syms[n]: v for n, v in provided_mags.items() if n in syms}
    eqs_subbed = [(orig, orig.subs(subs)) for orig in equations]

    # 3. Check consistency of fully-determined equations (no free symbols).
    # When all symbols in an equation are substituted, sympy simplifies
    # Eq(a, b) to the boolean atom sympy.true or sympy.false rather than
    # leaving it as an Eq object.
    conflicts: list[tuple[str, Any, Any]] = []
    eqs_with_unknowns: list[sympy.Eq] = []

    for orig_eq, subbed_eq in eqs_subbed:
        # Case 1: sympy already reduced to a boolean atom.
        # BooleanTrue means exact symbolic equality — no conflict.
        # BooleanFalse means exact symbolic inequality, but we still apply
        # a tolerance check to avoid false positives from floating-point
        # rounding (e.g. different paths to compute the same physical value
        # may differ by ~1 ULP).
        if subbed_eq is sympy.true:
            continue  # consistent, nothing to do
        if subbed_eq is sympy.false:
            lhs_sym = orig_eq.lhs
            lhs_name = (
                lhs_sym.name if isinstance(lhs_sym, sympy.Symbol) else str(lhs_sym)
            )
            try:
                lhs_val = float(sympy.N(orig_eq.lhs.subs(subs)))
                rhs_val = float(sympy.N(orig_eq.rhs.subs(subs)))
                scale = max(abs(rhs_val), 1e-100)
                if abs(lhs_val - rhs_val) / scale > tolerance:
                    provided_lhs = provided_mags.get(lhs_name, lhs_val)
                    conflicts.append((lhs_name, provided_lhs, rhs_val))
            except (TypeError, ValueError):
                pass  # can't evaluate numerically; skip
            continue

        # Case 2: still an Eq object — check free symbols
        free = subbed_eq.free_symbols
        if not free:
            # Fully numeric Eq that wasn't simplified to True/False
            try:
                lhs_val = float(sympy.N(subbed_eq.lhs))
                rhs_val = float(sympy.N(subbed_eq.rhs))
                scale = max(abs(rhs_val), 1e-100)
                if abs(lhs_val - rhs_val) / scale > tolerance:
                    lhs_sym = orig_eq.lhs
                    lhs_name = (
                        lhs_sym.name if isinstance(lhs_sym, sympy.Symbol) else str(lhs_sym)
                    )
                    provided_lhs = provided_mags.get(lhs_name, lhs_val)
                    conflicts.append((lhs_name, provided_lhs, rhs_val))
            except (TypeError, ValueError):
                eqs_with_unknowns.append(subbed_eq)
        else:
            eqs_with_unknowns.append(subbed_eq)

    if conflicts:
        raise ParameterConflict(conflicts)

    # 4. Solve for remaining unknowns (best-effort: don't raise if unsolvable)
    solved_mags = dict(provided_mags)

    if eqs_with_unknowns:
        all_unknowns = set()
        for eq in eqs_with_unknowns:
            all_unknowns.update(eq.free_symbols)
        unknowns_list = sorted(all_unknowns, key=lambda s: s.name)

        try:
            solutions = sympy.solve(eqs_with_unknowns, unknowns_list, dict=True)
        except Exception:
            solutions = []

        if solutions:
            solution = solutions[0]
            for sym, val in solution.items():
                try:
                    solved_mags[sym.name] = float(sympy.N(val))
                except (TypeError, ValueError):
                    pass

    return solved_mags


# ─── Public API ──────────────────────────────────────────────────────────────


def resolve_micro_params(
    provided: dict[str, Any],
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    """Resolve algebraic relationships among MicroParameters.

    Given a partial set of micro parameter values (which may include dependent
    parameters in place of independent ones), solves for any algebraically
    determinable values and returns an augmented dict.

    This function does **not** fill in defaults — call
    :func:`~lysis.config.paramcheck.load_micro_params` afterwards to supply
    defaults for any missing independent parameters and perform final
    consistency validation.

    The typical inversion case: a researcher specifies ``unbind_rate_tPA_woPLG``
    and ``bind_rate_tPA`` in the CSV; this function solves for
    ``diss_const_tPA_woPLG = unbind_rate_tPA_woPLG / bind_rate_tPA``.

    :param provided: Dict mapping parameter names to values.  Values may be
        :class:`~pint.Quantity` objects, unit-bearing strings (e.g. ``"72.7 nm"``),
        or plain numbers.  May include both independent and dependent
        MicroParameters fields.  Non-parameter keys are passed through unchanged.
    :type provided: dict[str, Quantity | str | int | float]
    :param tolerance: Relative tolerance for numerical consistency checks.
    :type tolerance: float
    :returns: Copy of *provided* augmented with any additionally solved
        parameter values, expressed as canonical-unit strings.
    :rtype: dict[str, str | int | float]
    :raises ParameterConflict: If fully-specified equations yield inconsistent
        values (e.g. ``bind_rate × diss_const ≠ unbind_rate``).
    """
    solved = _resolve_equations(
        _micro_equations(), _MICRO_SYMS, _MICRO_EQ_PARAMS, provided, tolerance
    )

    # Merge: start from original (preserves non-equation params), then overlay
    # solved values as canonical strings so load_micro_params can parse them.
    result = dict(provided)
    for name, magnitude in solved.items():
        result[name] = _mag_to_str(name, magnitude)

    return result


def resolve_macro_params(
    provided: dict[str, Any],
    micro_params_basedict: dict[str, Any],
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    """Resolve algebraic relationships among MacroParameters.

    Handles both the grid geometry equations and ``average_bound_time``, which
    depends on the microscale unbinding rate.

    :param provided: Dict of macro parameter values (same conventions as
        :func:`resolve_micro_params`).
    :type provided: dict[str, Quantity | str | int | float]
    :param micro_params_basedict: Fully resolved micro parameters as a base-type
        dict (e.g. from :meth:`~.parameters.MicroParameters.to_basedict`).
        Provides ``unbind_rate_tPA_woPLG`` for ``average_bound_time`` resolution.
    :type micro_params_basedict: dict[str, str | int | float]
    :param tolerance: Relative tolerance for numerical consistency checks.
    :type tolerance: float
    :returns: Copy of *provided* augmented with any additionally solved values.
    :rtype: dict[str, str | int | float]
    :raises ParameterConflict: If provided values are mutually inconsistent.
    """
    solved = _resolve_equations(
        _macro_equations(), _MACRO_SYMS, _MACRO_EQ_PARAMS, provided, tolerance
    )

    # Handle average_bound_time = 1 / unbind_rate_tPA_woPLG
    # (depends on a micro param, so it's not in the macro equation system)
    unbind_str = micro_params_basedict.get("unbind_rate_tPA_woPLG")
    if unbind_str is not None:
        try:
            unbind_mag = _parse_to_magnitude("unbind_rate_tPA_woPLG", unbind_str)
            computed_abt = 1.0 / unbind_mag  # canonical: seconds

            if "average_bound_time" in provided and provided["average_bound_time"] not in ("", None):
                provided_abt = _parse_to_magnitude(
                    "average_bound_time", provided["average_bound_time"]
                )
                scale = max(abs(computed_abt), 1e-100)
                if abs(provided_abt - computed_abt) / scale > tolerance:
                    raise ParameterConflict(
                        [("average_bound_time", provided_abt, computed_abt)]
                    )

            solved["average_bound_time"] = computed_abt
        except ParameterConflict:
            raise
        except Exception:
            pass  # If we can't compute it, leave it for MacroParameters.__post_init__

    result = dict(provided)
    for name, magnitude in solved.items():
        result[name] = _mag_to_str(name, magnitude)

    return result
