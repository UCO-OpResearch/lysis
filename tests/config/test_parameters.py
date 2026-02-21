"""Unit tests for lysis.config.parameters.

Tests cover:
    Parameters (base class)  -- units(), fortran_names(), to_dict(),
                                 to_basedict(), parse_from_basedict(),
                                 print_default_values(), __str__()
    MicroParameters          -- default/custom init, dependent parameter
                                 calculations, immutability, serialization
    MacroParameters          -- default/custom init, dependent parameter
                                 calculations, immutability, serialization

All physical assertions use pytest.approx() with default tolerances unless a
tighter or looser tolerance is explicitly justified.
"""

import dataclasses
import warnings

import pytest
from pint import Quantity

from lysis.config.constants import Q_
from lysis.config.parameters import MacroParameters, MicroParameters, Parameters


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _micro(**kwargs) -> MicroParameters:
    """Return a MicroParameters instance, optionally overriding fields."""
    return MicroParameters(**kwargs)


def _macro(micro=None, **kwargs) -> MacroParameters:
    """Return a MacroParameters instance.

    Suppresses RuntimeWarning that parse_from_basedict() emits when it
    encounters bare numeric values for unit-bearing parameters.
    """
    if micro is None:
        micro = _micro()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return MacroParameters(micro_params=micro, **kwargs)


# ===========================================================================
# MicroParameters — instantiation
# ===========================================================================


class TestMicroParametersDefaults:
    """MicroParameters instantiates correctly with all default values."""

    def test_instantiation_succeeds(self):
        """Default MicroParameters() raises no exception."""
        micro = _micro()
        assert isinstance(micro, MicroParameters)

    def test_inherits_parameters(self):
        """MicroParameters is a subclass of Parameters."""
        assert issubclass(MicroParameters, Parameters)


# ===========================================================================
# MicroParameters — custom initialisation
# ===========================================================================


class TestMicroParametersCustomInit:
    """Overriding independent parameters propagates correctly."""

    def test_override_fiber_radius(self):
        """fiber_radius can be overridden with a Quantity."""
        micro = _micro(fiber_radius=Q_("50 nm"))
        assert micro.fiber_radius.to("nm").magnitude == pytest.approx(50.0)

    def test_override_nodes_in_micro_row(self):
        """nodes_in_micro_row can be overridden with an int."""
        micro = _micro(nodes_in_micro_row=13)
        assert micro.nodes_in_micro_row == 13

    def test_override_micro_simulations(self):
        """micro_simulations can be overridden."""
        micro = _micro(micro_simulations=1_000)
        assert micro.micro_simulations == 1_000

    def test_override_micro_seed(self):
        """micro_seed can be overridden."""
        micro = _micro(micro_seed=42)
        assert micro.micro_seed == 42

    def test_override_diss_const_tPA_woPLG(self):
        """diss_const_tPA_woPLG override propagates into unbind_rate_tPA_woPLG."""
        new_kd = Q_("1.0 micromolar")
        micro = _micro(
            diss_const_tPA_woPLG=new_kd, bind_rate_tPA=Q_("0.1 (micromolar*sec)^-1")
        )
        assert micro.unbind_rate_tPA_woPLG.to("sec^-1").magnitude == pytest.approx(0.1)


# ===========================================================================
# MicroParameters — immutability
# ===========================================================================


class TestMicroParametersImmutability:
    """MicroParameters is frozen; attribute assignment must raise FrozenInstanceError."""

    def test_cannot_set_independent_param(self):
        """Setting an independent param on a frozen instance raises."""
        micro = _micro()
        with pytest.raises(Exception):  # dataclasses.FrozenInstanceError
            micro.nodes_in_micro_row = 99

    def test_cannot_set_dependent_param(self):
        """Setting a dependent param on a frozen instance raises."""
        micro = _micro()
        with pytest.raises(Exception):
            micro.protofibril_radius = Q_("1 nm")


# ===========================================================================
# MicroParameters — serialisation
# ===========================================================================


class TestMicroParametersSerialization:
    """to_dict(), to_basedict(), __str__(), and parse_from_basedict()."""

    def test_to_dict_returns_dict(self):
        """to_dict() returns a dict."""
        assert isinstance(_micro().to_dict(), dict)

    def test_to_dict_contains_fiber_radius_as_quantity(self):
        """to_dict() preserves Quantity types (not strings)."""
        d = _micro().to_dict()
        assert isinstance(d["fiber_radius"], Quantity)

    def test_to_dict_contains_nodes_in_micro_row(self):
        """to_dict() includes the int-valued independent parameter."""
        d = _micro().to_dict()
        assert d["nodes_in_micro_row"] == 7

    def test_to_basedict_returns_dict(self):
        """to_basedict() returns a dict."""
        assert isinstance(_micro().to_basedict(), dict)

    def test_to_basedict_converts_quantities_to_strings(self):
        """to_basedict() converts Quantity values to strings."""
        bd = _micro().to_basedict()
        assert isinstance(bd["fiber_radius"], str)
        assert isinstance(bd["diss_const_tPA_woPLG"], str)

    def test_to_basedict_passes_through_ints(self):
        """to_basedict() passes through int values unchanged."""
        bd = _micro().to_basedict()
        assert bd["nodes_in_micro_row"] == 7
        assert isinstance(bd["nodes_in_micro_row"], int)

    def test_to_basedict_passes_through_floats(self):
        """to_basedict() passes through float values unchanged."""
        bd = _micro().to_basedict()
        assert isinstance(bd["snap_proportion"], float)

    def test_to_basedict_quantity_string_parseable(self):
        """Quantity strings in to_basedict() can be parsed back by Pint."""
        bd = _micro().to_basedict()
        reparsed = Q_(bd["fiber_radius"])
        assert reparsed.to("nm").magnitude == pytest.approx(72.7 / 2)

    def test_str_returns_string(self):
        """__str__() returns a non-empty string."""
        s = str(_micro())
        assert isinstance(s, str)
        assert len(s) > 0

    def test_print_default_values_returns_string(self):
        """print_default_values() returns a non-empty string."""
        s = MicroParameters.print_default_values()
        assert isinstance(s, str)
        assert len(s) > 0

    def test_parse_from_basedict_roundtrip(self, capsys):
        """parse_from_basedict(to_basedict()) reproduces independent params."""
        micro = _micro()
        bd = micro.to_basedict()
        loaded = MicroParameters.parse_from_basedict(bd)
        capsys.readouterr()  # suppress debug print

        assert loaded.nodes_in_micro_row == micro.nodes_in_micro_row
        assert loaded.micro_simulations == micro.micro_simulations
        assert loaded.micro_seed == micro.micro_seed
        assert loaded.fiber_radius.to("nm").magnitude == pytest.approx(
            micro.fiber_radius.to("nm").magnitude
        )

    def test_parse_from_basedict_ignores_dependent_params(self, capsys):
        """parse_from_basedict() ignores stored dependent parameters silently."""
        micro = _micro()
        bd = micro.to_basedict()
        # Tamper a dependent param; it should be ignored and recalculated
        bd["protofibril_radius"] = "999 nm"
        loaded = MicroParameters.parse_from_basedict(bd)
        capsys.readouterr()
        # Recalculated value should match the original, not the tampered value
        assert loaded.protofibril_radius.to("nm").magnitude == pytest.approx(
            micro.protofibril_radius.to("nm").magnitude
        )

    def test_parse_from_basedict_ignores_unknown_keys(self, capsys):
        """parse_from_basedict() silently ignores keys absent from the constructor."""
        micro = _micro()
        bd = micro.to_basedict()
        bd["nonexistent_key"] = "some_value"
        loaded = MicroParameters.parse_from_basedict(bd)
        capsys.readouterr()
        assert loaded.nodes_in_micro_row == micro.nodes_in_micro_row

    def test_parse_from_basedict_warns_on_bare_number_for_units_param(self, capsys):
        """parse_from_basedict() issues RuntimeWarning for bare numeric unit-bearing params."""
        micro = _micro()
        bd = micro.to_basedict()
        # Replace the string Quantity with a raw float
        bd["fiber_radius"] = 36.35  # bare number, no units
        with pytest.warns(RuntimeWarning, match="fiber_radius"):
            MicroParameters.parse_from_basedict(bd)
        capsys.readouterr()


# ===========================================================================
# MicroParameters — metadata (units, fortran_names)
# ===========================================================================


class TestMicroParametersMetadata:
    """units() and fortran_names() parse docstrings correctly."""

    def test_units_returns_dict(self):
        """units() returns a dict."""
        assert isinstance(MicroParameters.units(), dict)

    def test_units_contains_fiber_radius(self):
        """fiber_radius is mapped to 'microns' in units()."""
        u = MicroParameters.units()
        assert "fiber_radius" in u
        assert u["fiber_radius"] == "microns"

    def test_units_excludes_none_params(self):
        """Parameters with ':Units: None' are not in the units() dict."""
        u = MicroParameters.units()
        assert "nodes_in_micro_row" not in u
        assert "micro_seed" not in u
        assert "snap_proportion" not in u

    def test_units_contains_rate_params(self):
        """Binding and unbinding rate parameters appear in units()."""
        u = MicroParameters.units()
        assert "bind_rate_tPA" in u
        assert "bind_rate_PLG" in u
        assert "unbind_rate_PLG_intact" in u
        assert "unbind_rate_tPA_woPLG" in u

    def test_fortran_names_returns_dict(self):
        """fortran_names() returns a dict."""
        assert isinstance(MicroParameters.fortran_names(), dict)

    def test_fortran_names_fiber_radius(self):
        """fiber_radius maps to 'radius' in Fortran."""
        fn = MicroParameters.fortran_names()
        assert fn.get("fiber_radius") == "radius"

    def test_fortran_names_diss_consts(self):
        """Dissociation constants have the expected Fortran names."""
        fn = MicroParameters.fortran_names()
        assert fn.get("diss_const_tPA_woPLG") == "KdtPAnoplg"
        assert fn.get("diss_const_tPA_wPLG") == "KdtPAyesplg"
        assert fn.get("diss_const_PLG_intact") == "KdPLGintact"
        assert fn.get("diss_const_PLG_nicked") == "KdPLGnicked"

    def test_fortran_names_excludes_none_params(self):
        """Parameters with ':Fortran: None' are not in fortran_names()."""
        fn = MicroParameters.fortran_names()
        assert "fibrinogen_length" not in fn
        assert "fibrinogen_radius" not in fn
        assert "protein_per_fiber" not in fn

    def test_fortran_names_micro_simulations(self):
        """micro_simulations maps to 'simulations' in Fortran."""
        fn = MicroParameters.fortran_names()
        assert fn.get("micro_simulations") == "simulations"

    def test_fortran_names_binding_sites(self):
        """binding_sites maps to 'bs' in Fortran."""
        fn = MicroParameters.fortran_names()
        assert fn.get("binding_sites") == "bs"


# ===========================================================================
# MacroParameters — instantiation
# ===========================================================================


class TestMacroParametersDefaults:
    """MacroParameters instantiates correctly with default values when given a MicroParameters."""

    def test_instantiation_succeeds(self):
        """MacroParameters(micro_params=...) raises no exception."""
        assert isinstance(_macro(), MacroParameters)

    def test_inherits_parameters(self):
        """MacroParameters is a subclass of Parameters."""
        assert issubclass(MacroParameters, Parameters)

    def test_requires_micro_params(self):
        """MacroParameters() without micro_params raises TypeError."""
        with pytest.raises(TypeError):
            MacroParameters()  # micro_params has no default

    def test_micro_params_stored(self):
        """micro_params attribute is the MicroParameters instance passed in."""
        micro = _micro()
        macro = _macro(micro)
        assert macro.micro_params is micro


# ===========================================================================
# MacroParameters — custom initialisation
# ===========================================================================


class TestMacroParametersCustomInit:
    """Overriding independent parameters propagates correctly."""

    def test_override_cols(self):
        """Overriding cols changes cols and all derived grid quantities."""
        macro = _macro(cols=50)
        assert macro.cols == 50
        assert macro.full_row == 3 * 50 - 1
        assert macro.xz_row == 2 * 50 - 1

    def test_override_rows_and_empty_rows(self):
        """Overriding rows and empty_rows changes fiber_rows."""
        macro = _macro(rows=100, empty_rows=20)
        assert macro.fiber_rows == 80

    def test_override_macro_seed(self):
        """Overriding macro_seed propagates to the RNG state tuple."""
        macro = _macro(macro_seed=12345)
        assert macro.state[3] == 12345

    def test_override_total_time(self):
        """Overriding total_time changes total_time_steps and number_of_saves."""
        # Use a simple total_time divisible by save_interval for exact arithmetic
        macro_a = _macro(total_time=Q_("100 sec"), save_interval=Q_("10 sec"))
        macro_b = _macro(total_time=Q_("200 sec"), save_interval=Q_("10 sec"))
        assert macro_b.number_of_saves > macro_a.number_of_saves

    def test_override_pore_size_changes_time_step(self):
        """Overriding pore_size changes time_step (Equation 2.4)."""
        macro_small = _macro(pore_size=Q_("0.5 um"))
        macro_large = _macro(pore_size=Q_("2.0 um"))
        # Larger pore → larger time step
        assert macro_large.time_step.magnitude > macro_small.time_step.magnitude


# ===========================================================================
# MacroParameters — immutability
# ===========================================================================


class TestMacroParametersImmutability:
    """MacroParameters is frozen; attribute assignment must raise."""

    def test_cannot_set_independent_param(self):
        """Setting an independent param on a frozen macro instance raises."""
        macro = _macro()
        with pytest.raises(Exception):
            macro.cols = 99

    def test_cannot_set_dependent_param(self):
        """Setting a dependent param on a frozen macro instance raises."""
        macro = _macro()
        with pytest.raises(Exception):
            macro.full_row = 999


# ===========================================================================
# MacroParameters — serialisation
# ===========================================================================


class TestMacroParametersSerialization:
    """to_dict(), to_basedict(), and parse_from_basedict() for MacroParameters."""

    def test_to_dict_returns_dict(self):
        """to_dict() returns a dict."""
        assert isinstance(_macro().to_dict(), dict)

    def test_to_basedict_returns_dict(self):
        """to_basedict() returns a dict."""
        assert isinstance(_macro().to_basedict(), dict)

    def test_to_basedict_skips_micro_params_dict(self):
        """to_basedict() skips the nested micro_params dict (it's a sub-dict)."""
        bd = _macro().to_basedict()
        # micro_params is a nested dict, so it should be skipped
        assert "micro_params" not in bd

    def test_to_basedict_converts_pore_size_to_string(self):
        """to_basedict() converts pore_size Quantity to a string."""
        bd = _macro().to_basedict()
        assert isinstance(bd["pore_size"], str)

    def test_to_basedict_passes_through_cols(self):
        """to_basedict() passes through int-valued cols unchanged."""
        bd = _macro(cols=93).to_basedict()
        assert bd["cols"] == 93

    def test_str_returns_non_empty_string(self):
        """__str__() returns a non-empty string for MacroParameters."""
        s = str(_macro())
        assert isinstance(s, str)
        assert len(s) > 0

    def test_print_default_values_raises_without_micro(self):
        """print_default_values() raises TypeError because micro_params has no default."""
        with pytest.raises(TypeError):
            MacroParameters.print_default_values()

    def test_parse_from_basedict_raises_without_micro_params(self, capsys):
        """parse_from_basedict() raises TypeError for MacroParameters.

        to_basedict() skips the micro_params nested dict, so it is absent
        from the basedict. parse_from_basedict() does not accept micro_params
        as a separate argument, so MacroParameters construction fails with
        TypeError because the required micro_params field is missing.
        """
        micro = _micro()
        macro = _macro(micro, cols=50, rows=100, empty_rows=10)
        bd = macro.to_basedict()
        # micro_params is excluded from to_basedict(); the roundtrip therefore
        # cannot succeed without an external mechanism (see paramcheck.load_macro_params).
        with pytest.raises(TypeError, match="micro_params"):
            MacroParameters.parse_from_basedict(bd)
        capsys.readouterr()


# ===========================================================================
# MicroParameters vs MacroParameters — key separation
# ===========================================================================


class TestParametersKeySeparation:
    """MicroParameters and MacroParameters must not share any field names."""

    def test_no_overlapping_field_names(self):
        """MicroParameters and MacroParameters define no fields with the same name."""
        micro_keys = {f.name for f in dataclasses.fields(MicroParameters)}
        macro_keys = {f.name for f in dataclasses.fields(MacroParameters)}
        overlap = micro_keys & macro_keys
        assert overlap == set(), f"Overlapping field names: {overlap}"
