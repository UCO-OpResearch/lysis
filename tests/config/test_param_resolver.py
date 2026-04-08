"""Unit tests for lysis.config.param_resolver.

Tests cover:
  - Forward computation: provide independent params, verify dependent params solved
  - Inversion: provide dependent param + some independents, solve for missing independent
  - Conflict detection: inconsistent over-constrained values
  - Pass-through: non-equation params are preserved unchanged
  - Unit handling: non-canonical units are converted transparently
  - Macro resolution: grid geometry and time_step equations
  - average_bound_time: depends on micro unbind rate
"""

import math

import pytest
from pint import Quantity

from lysis.config.constants import Q_
from lysis.config.param_resolver import (
    ParameterConflict,
    UnderdeterminedParameters,
    resolve_macro_params,
    resolve_micro_params,
)
from lysis.config.parameters import MacroParameters, MicroParameters


# ─── Helpers ─────────────────────────────────────────────────────────────────


def _micro_defaults() -> dict:
    """Return to_basedict() for default MicroParameters."""
    return MicroParameters().to_basedict()


def _macro_defaults(micro=None) -> dict:
    """Return to_basedict() for default MacroParameters."""
    mp = micro or MicroParameters()
    return MacroParameters(micro_params=mp).to_basedict()


def _parse_q(s) -> float:
    """Parse a canonical string to float magnitude."""
    return float(Q_(s).magnitude) if isinstance(s, str) else float(s)


# ─── resolve_micro_params: forward pass ──────────────────────────────────────


class TestResolveмикроForward:
    """Providing all independent params → dependent params computed correctly."""

    def test_protofibril_radius_computed(self):
        defaults = _micro_defaults()
        result = resolve_micro_params(defaults)
        expected = float(Q_("2.4 nanometers").to("microns").magnitude)
        got = _parse_q(result["protofibril_radius"])
        assert math.isclose(got, expected, rel_tol=1e-6)

    def test_unbind_rate_PLG_intact_computed(self):
        defaults = _micro_defaults()
        result = resolve_micro_params(defaults)
        # koff = kon × Kd = 0.1 × 38 = 3.8 sec^-1
        got = _parse_q(result["unbind_rate_PLG_intact"])
        assert math.isclose(got, 3.8, rel_tol=1e-6)

    def test_unbind_rate_tPA_woPLG_computed(self):
        defaults = _micro_defaults()
        result = resolve_micro_params(defaults)
        # 0.1 × 0.36 = 0.036 sec^-1
        got = _parse_q(result["unbind_rate_tPA_woPLG"])
        assert math.isclose(got, 0.036, rel_tol=1e-6)

    def test_fibrin_conc_per_fiber_computed(self):
        defaults = _micro_defaults()
        result = resolve_micro_params(defaults)
        # Default ~871 µM from Bannish et al. 2017
        got = _parse_q(result["fibrin_conc_per_fiber"])
        assert 800 < got < 1000, f"Expected ~871 µM, got {got}"

    def test_binding_sites_computed(self):
        defaults = _micro_defaults()
        result = resolve_micro_params(defaults)
        # Default ~427 µM from Bannish et al. 2017
        got = _parse_q(result["binding_sites"])
        assert 400 < got < 500, f"Expected ~427 µM, got {got}"

    def test_non_equation_params_pass_through(self):
        defaults = _micro_defaults()
        defaults["micro_simulations"] = 12345
        defaults["micro_version"] = "test_version"
        result = resolve_micro_params(defaults)
        assert result["micro_simulations"] == 12345
        assert result["micro_version"] == "test_version"


# ─── resolve_micro_params: inversion ─────────────────────────────────────────


class TestResolveMicroInversion:
    """Providing a dependent param + one independent → missing independent solved."""

    def test_solve_diss_const_from_unbind_rate(self):
        """Provide unbind_rate_tPA_woPLG + bind_rate_tPA → solve diss_const_tPA_woPLG."""
        provided = {
            "bind_rate_tPA": "0.1 1 / micromolar / second",
            "unbind_rate_tPA_woPLG": "0.072 1 / second",
            # diss_const_tPA_woPLG is NOT provided; should be solved as 0.72 µM
        }
        result = resolve_micro_params(provided)
        assert "diss_const_tPA_woPLG" in result
        got = _parse_q(result["diss_const_tPA_woPLG"])
        assert math.isclose(got, 0.72, rel_tol=1e-6), f"Expected 0.72, got {got}"

    def test_solve_bind_rate_from_unbind_rate(self):
        """Provide unbind_rate_PLG_intact + diss_const → solve bind_rate_PLG."""
        provided = {
            "diss_const_PLG_intact": "38 micromolar",
            "unbind_rate_PLG_intact": "3.8 1 / second",
        }
        result = resolve_micro_params(provided)
        assert "bind_rate_PLG" in result
        got = _parse_q(result["bind_rate_PLG"])
        assert math.isclose(got, 0.1, rel_tol=1e-6)

    def test_solve_fibrinogen_radius_from_protofibril(self):
        """protofibril_radius = 2 × fibrinogen_radius → solve fibrinogen_radius."""
        provided = {"protofibril_radius": "3.0 nanometers"}
        result = resolve_micro_params(provided)
        assert "fibrinogen_radius" in result
        got = _parse_q(result["fibrinogen_radius"])
        # 3.0 nm / 2 = 1.5 nm; in microns:
        expected = float(Q_("1.5 nanometers").to("microns").magnitude)
        assert math.isclose(got, expected, rel_tol=1e-6)

    def test_non_canonical_units_accepted(self):
        """Values in non-canonical units (nm vs µm) are handled correctly."""
        # fiber_radius canonical unit is microns; provide in nm
        provided = {"fiber_radius": "72.7 nm"}
        result = resolve_micro_params(provided)
        # Should be stored in canonical units (microns)
        got = _parse_q(result["fiber_radius"])
        expected = float(Q_("72.7 nm").to("microns").magnitude)
        assert math.isclose(got, expected, rel_tol=1e-6)

    def test_pint_quantity_input_accepted(self):
        """Pint Quantity objects (not strings) are accepted."""
        provided = {"fiber_radius": Q_("72.7 nm")}
        result = resolve_micro_params(provided)
        got = _parse_q(result["fiber_radius"])
        expected = float(Q_("72.7 nm").to("microns").magnitude)
        assert math.isclose(got, expected, rel_tol=1e-6)


# ─── resolve_micro_params: conflict detection ────────────────────────────────


class TestResolveMicroConflict:
    """Inconsistent over-constrained values → ParameterConflict raised."""

    def test_conflict_koff_ne_kon_times_kd(self):
        """bind_rate × diss_const ≠ unbind_rate → ParameterConflict."""
        provided = {
            "bind_rate_tPA": "0.1 1 / micromolar / second",
            "diss_const_tPA_woPLG": "0.36 micromolar",
            "unbind_rate_tPA_woPLG": "9.99 1 / second",  # inconsistent: should be 0.036
        }
        with pytest.raises(ParameterConflict) as exc_info:
            resolve_micro_params(provided)
        assert any(
            "unbind_rate_tPA_woPLG" in c[0] for c in exc_info.value.conflicts
        )

    def test_conflict_protofibril_ne_2_fibrinogen(self):
        """protofibril_radius ≠ 2 × fibrinogen_radius → ParameterConflict."""
        provided = {
            "fibrinogen_radius": "1.2 nanometers",
            "protofibril_radius": "5.0 nanometers",  # should be 2.4 nm
        }
        with pytest.raises(ParameterConflict) as exc_info:
            resolve_micro_params(provided)
        assert any(
            "protofibril_radius" in c[0] for c in exc_info.value.conflicts
        )

    def test_conflict_reports_all_at_once(self):
        """Multiple conflicts are all reported, not just the first one."""
        provided = {
            "bind_rate_tPA": "0.1 1 / micromolar / second",
            "diss_const_tPA_woPLG": "0.36 micromolar",
            "unbind_rate_tPA_woPLG": "9.99 1 / second",  # wrong
            "fibrinogen_radius": "1.2 nanometers",
            "protofibril_radius": "5.0 nanometers",  # wrong
        }
        with pytest.raises(ParameterConflict) as exc_info:
            resolve_micro_params(provided)
        # At least two conflicts should be reported
        assert len(exc_info.value.conflicts) >= 2


# ─── resolve_macro_params: forward pass ──────────────────────────────────────


class TestResolveMacroForward:
    """Providing all independent macro params → dependent params computed correctly."""

    def setup_method(self):
        self.micro = MicroParameters()
        self.micro_dict = self.micro.to_basedict()
        self.macro_defaults = _macro_defaults(self.micro)

    def test_full_row_computed(self):
        result = resolve_macro_params(self.macro_defaults, self.micro_dict)
        expected = 3 * 93 - 1  # cols=93 → full_row=278
        got = result["full_row"]
        assert got == expected

    def test_xz_row_computed(self):
        result = resolve_macro_params(self.macro_defaults, self.micro_dict)
        expected = 2 * 93 - 1  # 185
        got = result["xz_row"]
        assert got == expected

    def test_total_edges_computed(self):
        result = resolve_macro_params(self.macro_defaults, self.micro_dict)
        expected = MacroParameters(micro_params=self.micro).total_edges
        got = result["total_edges"]
        assert got == expected

    def test_time_step_computed(self):
        result = resolve_macro_params(self.macro_defaults, self.micro_dict)
        expected = float(MacroParameters(micro_params=self.micro).time_step.to("seconds").magnitude)
        got = _parse_q(result["time_step"])
        assert math.isclose(got, expected, rel_tol=1e-6)

    def test_average_bound_time_computed(self):
        result = resolve_macro_params(self.macro_defaults, self.micro_dict)
        # 1 / unbind_rate_tPA_woPLG = 1 / 0.036 ≈ 27.78 sec
        got = _parse_q(result["average_bound_time"])
        expected = float(
            (1.0 / Q_("0.036 1/second")).to("seconds").magnitude
        )
        assert math.isclose(got, expected, rel_tol=1e-4)


# ─── resolve_macro_params: inversion ────────────────────────────────────────


class TestResolveMacroInversion:
    """Providing derived geometry params → solve for missing independent."""

    def setup_method(self):
        self.micro = MicroParameters()
        self.micro_dict = self.micro.to_basedict()

    def test_solve_cols_from_full_row(self):
        """full_row = 3*cols - 1 → solve for cols."""
        provided = {"full_row": 278}  # → cols = (278 + 1) / 3 = 93
        result = resolve_macro_params(provided, self.micro_dict)
        assert "cols" in result
        assert result["cols"] == 93

    def test_solve_empty_rows_from_fiber_rows(self):
        """fiber_rows = rows - empty_rows → solve for empty_rows."""
        provided = {"rows": 121, "fiber_rows": 93}  # → empty_rows = 28
        result = resolve_macro_params(provided, self.micro_dict)
        assert "empty_rows" in result
        assert result["empty_rows"] == 28

    def test_solve_moving_probability_from_time_step(self):
        """time_step = q * dx^2 / (12 * D) → solve for moving_probability."""
        # Default time_step from MacroParameters
        macro = MacroParameters(micro_params=self.micro)
        dt = float(macro.time_step.to("seconds").magnitude)
        dx = float(Q_("1.0135 um").to("centimeters").magnitude)
        D = 5.0e-7  # cm^2/s (canonical)
        provided = {
            "time_step": f"{dt} seconds",
            "pore_size": "1.0135e-4 centimeters",
            "diffusion_coeff": "5e-7 cm^2/s",
        }
        result = resolve_macro_params(provided, self.micro_dict)
        assert "moving_probability" in result
        got = _parse_q(result["moving_probability"])
        assert math.isclose(got, 0.2, rel_tol=1e-4)


# ─── resolve_macro_params: conflict detection ────────────────────────────────


class TestResolveMacroConflict:
    """Inconsistent macro values → ParameterConflict raised."""

    def setup_method(self):
        self.micro = MicroParameters()
        self.micro_dict = self.micro.to_basedict()

    def test_conflict_full_row_inconsistent(self):
        """full_row ≠ 3*cols - 1 → ParameterConflict."""
        provided = {"cols": 93, "full_row": 999}  # should be 278
        with pytest.raises(ParameterConflict) as exc_info:
            resolve_macro_params(provided, self.micro_dict)
        assert any("full_row" in c[0] for c in exc_info.value.conflicts)

    def test_conflict_average_bound_time_inconsistent(self):
        """average_bound_time ≠ 1/unbind_rate → ParameterConflict."""
        provided = {"average_bound_time": "999 seconds"}
        with pytest.raises(ParameterConflict) as exc_info:
            resolve_macro_params(provided, self.micro_dict)
        assert any("average_bound_time" in c[0] for c in exc_info.value.conflicts)


# ─── Round-trip integration ──────────────────────────────────────────────────


class TestRoundTrip:
    """Resolved dicts can be passed to load_micro_params / load_macro_params."""

    def test_micro_round_trip(self):
        from lysis.config.paramcheck import load_micro_params

        defaults = _micro_defaults()
        resolved = resolve_micro_params(defaults)
        # load_micro_params should succeed without error
        loaded = load_micro_params(resolved)
        assert isinstance(loaded, MicroParameters)
        # Spot-check a value
        assert math.isclose(
            loaded.fiber_radius.to("microns").magnitude,
            MicroParameters().fiber_radius.to("microns").magnitude,
            rel_tol=1e-6,
        )

    def test_macro_round_trip(self):
        from lysis.config.paramcheck import load_macro_params, load_micro_params

        micro = MicroParameters()
        micro_dict = micro.to_basedict()
        macro_dict = _macro_defaults(micro)

        resolved_micro = resolve_micro_params(micro_dict)
        micro_loaded = load_micro_params(resolved_micro)

        resolved_macro = resolve_macro_params(macro_dict, micro_loaded.to_basedict())
        macro_loaded = load_macro_params(resolved_macro, micro_loaded)
        assert isinstance(macro_loaded, MacroParameters)
        assert macro_loaded.cols == 93
        assert macro_loaded.rows == 121

    def test_micro_inversion_round_trip(self):
        """Solving bind_rate from unbind_rate produces a valid MicroParameters."""
        from lysis.config.paramcheck import load_micro_params

        # Provide unbind_rate_tPA_woPLG instead of bind_rate_tPA
        provided = {
            "diss_const_tPA_woPLG": "0.72 micromolar",
            "unbind_rate_tPA_woPLG": "0.072 1 / second",
            # bind_rate_tPA is NOT provided → should be solved as 0.1
        }
        resolved = resolve_micro_params(provided)
        # bind_rate_tPA should be solved
        assert "bind_rate_tPA" in resolved
        got = _parse_q(resolved["bind_rate_tPA"])
        assert math.isclose(got, 0.1, rel_tol=1e-5)

        # load_micro_params should succeed (fills remaining defaults)
        loaded = load_micro_params(resolved)
        assert isinstance(loaded, MicroParameters)
        assert math.isclose(
            loaded.bind_rate_tPA.to("1/(micromolar*second)").magnitude,
            0.1,
            rel_tol=1e-5,
        )
