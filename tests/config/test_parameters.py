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

import numpy as np
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

    def test_micro_seed_default_is_uint32(self):
        """Default micro_seed has dtype np.uint32."""
        micro = _micro()
        assert isinstance(micro.micro_seed, np.uint32)

    def test_micro_seed_override_normalised_to_uint32(self):
        """Plain-int override is normalised to np.uint32."""
        micro = _micro(micro_seed=42)
        assert isinstance(micro.micro_seed, np.uint32)
        assert micro.micro_seed == np.uint32(42)

    def test_micro_seed_negative_wraps_bitwise(self):
        """Negative int override is reinterpreted as the uint32 bit pattern."""
        micro = _micro(micro_seed=-1)
        assert isinstance(micro.micro_seed, np.uint32)
        assert micro.micro_seed == np.uint32(0xFFFFFFFF)

    def test_micro_seed_high_bit_preserved(self):
        """High-bit int override is preserved bit-for-bit as uint32."""
        micro = _micro(micro_seed=0x80000001)
        assert isinstance(micro.micro_seed, np.uint32)
        assert micro.micro_seed == np.uint32(0x80000001)

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

    def test_fortran_names_micro_seed_has_uint32_suffix(self):
        """micro_seed carries the `|uint32` tag-direction suffix."""
        fn = MicroParameters.fortran_names()
        assert fn.get("micro_seed") == "seed|uint32"

    def test_inverse_fortran_map_decodes_uint32_suffix(self):
        """inverse_fortran_map() strips the suffix and records the transform."""
        inverse = MicroParameters.inverse_fortran_map()
        assert "seed" in inverse
        base_name, transform, cls = inverse["seed"]
        assert transform == "uint32"
        assert base_name == "micro_seed"

    def test_apply_fortran_transform_uint32(self):
        """apply_fortran_transform() reinterprets Fortran int32 as np.uint32 bits."""
        # Fortran's signed print of 0xDEADBEEF is -559038737.
        result = Parameters.apply_fortran_transform(-559038737, "uint32")
        assert isinstance(result, np.uint32)
        assert result == np.uint32(0xDEADBEEF)


class TestMetadataRegexCompleteness:
    """Verify units() and fortran_names() capture every tagged parameter.

    These tests independently parse parameters.py to find all :Units: and
    :Fortran: tags, then check that the regex-based methods don't silently
    drop any of them (e.g. because a docstring contains double-quotes).
    """

    @staticmethod
    def _parse_tags(tag_name):
        """Parse all non-None :tag_name: values from parameter docstrings."""
        import pkgutil
        import re

        text = pkgutil.get_data("lysis.config", "parameters.py").decode("utf-8")
        # Match 4-space-indented field declarations whose docstring contains
        # the given tag.  The lookahead prevents crossing """ boundaries.
        tag_pattern = re.compile(
            r"^\s{4}([a-zA-Z0-9_]+):.*\n"  # field declaration
            r"\s{4}\"{3}(?:(?!\"{3})[\s\S])*?"  # docstring content (stops at """)
            rf":{tag_name}:\s+([^\s\"]+)",  # tag value (excludes quotes)
            re.M,
        )
        return {
            m.group(1): m.group(2)
            for m in tag_pattern.finditer(text)
            if m.group(2) != "None"
        }

    def test_no_missing_fortran_tags(self):
        """Every non-None :Fortran: tag in parameters.py is captured."""
        expected = self._parse_tags("Fortran")
        actual = MicroParameters.fortran_names()
        missing = set(expected) - set(actual)
        assert not missing, (
            f"fortran_names() is missing parameters: {missing}. "
            f"The regex likely fails on docstrings containing double-quotes."
        )

    def test_no_missing_units_tags(self):
        """Every non-None :Units: tag in parameters.py is captured."""
        expected = self._parse_tags("Units")
        actual = MicroParameters.units()
        missing = set(expected) - set(actual)
        assert not missing, (
            f"units() is missing parameters: {missing}. "
            f"The regex likely fails on docstrings containing double-quotes."
        )


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

    def test_macro_seed_default_is_uint32(self):
        """Default macro_seed has dtype np.uint32."""
        macro = _macro()
        assert isinstance(macro.macro_seed, np.uint32)

    def test_macro_seed_negative_wraps_bitwise(self):
        """Negative int override is reinterpreted as the uint32 bit pattern."""
        macro = _macro(macro_seed=-1)
        assert isinstance(macro.macro_seed, np.uint32)
        assert macro.macro_seed == np.uint32(0xFFFFFFFF)
        # The RNG state tuple mirrors the normalised value.
        assert macro.state[3] == np.uint32(0xFFFFFFFF)

    def test_macro_seed_high_bit_preserved(self):
        """High-bit int override is preserved bit-for-bit as uint32."""
        macro = _macro(macro_seed=0xDEADBEEF)
        assert isinstance(macro.macro_seed, np.uint32)
        assert macro.macro_seed == np.uint32(0xDEADBEEF)

    def test_macro_state_tuple_elements_are_uint32(self):
        """All four RNG state entries are np.uint32."""
        macro = _macro(macro_seed=0xDEADBEEF)
        for element in macro.state:
            assert isinstance(element, np.uint32)

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
# MacroParameters — forced_unbind default and calculate_forced_unbind()
# ===========================================================================


class TestMacroParametersForcedUnbind:
    """forced_unbind field and MacroParameters.calculate_forced_unbind()."""

    def test_forced_unbind_default_is_nan(self):
        """Default forced_unbind is float('nan') to prevent accidental use."""
        import math

        macro = _macro()
        assert math.isnan(macro.forced_unbind)

    def test_calculate_forced_unbind_basic(self):
        """calculate_forced_unbind returns the correct fraction."""
        import numpy as np

        pli = np.array([True, True, False, False])
        kin = np.array([False, False, True, False])
        # 2 forced, 1 kinetic → 2/3
        result = MacroParameters.calculate_forced_unbind(pli, kin)
        assert result == pytest.approx(2 / 3)

    def test_calculate_forced_unbind_all_forced(self):
        """All events forced → result is 1.0."""
        import numpy as np

        pli = np.ones(10, dtype=bool)
        kin = np.zeros(10, dtype=bool)
        assert MacroParameters.calculate_forced_unbind(pli, kin) == pytest.approx(1.0)

    def test_calculate_forced_unbind_all_kinetic(self):
        """All events kinetic → result is 0.0."""
        import numpy as np

        pli = np.zeros(10, dtype=bool)
        kin = np.ones(10, dtype=bool)
        assert MacroParameters.calculate_forced_unbind(pli, kin) == pytest.approx(0.0)

    def test_calculate_forced_unbind_no_events_raises(self):
        """Both arrays all-False → raises ValueError."""
        import numpy as np

        pli = np.zeros(10, dtype=bool)
        kin = np.zeros(10, dtype=bool)
        with pytest.raises(ValueError, match="no unbinding events"):
            MacroParameters.calculate_forced_unbind(pli, kin)


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


# ===========================================================================
# Parameters.class_fortran_names()
# ===========================================================================


class TestClassFortranNames:
    """class_fortran_names() filters fortran_names() to the calling class."""

    def test_micro_only_contains_micro_fields(self):
        """MicroParameters.class_fortran_names() only has MicroParameters fields."""
        cfn = MicroParameters.class_fortran_names()
        micro_fields = {f.name for f in dataclasses.fields(MicroParameters)}
        assert set(cfn.keys()) <= micro_fields

    def test_macro_only_contains_macro_fields(self):
        """MacroParameters.class_fortran_names() only has MacroParameters fields."""
        cfn = MacroParameters.class_fortran_names()
        macro_fields = {f.name for f in dataclasses.fields(MacroParameters)}
        assert set(cfn.keys()) <= macro_fields

    def test_micro_contains_fiber_radius(self):
        """fiber_radius is a micro field and should be present."""
        cfn = MicroParameters.class_fortran_names()
        assert "fiber_radius" in cfn
        assert cfn["fiber_radius"] == "radius"

    def test_macro_contains_cols(self):
        """cols is a macro field and should be present."""
        cfn = MacroParameters.class_fortran_names()
        assert "cols" in cfn

    def test_micro_excludes_macro_fields(self):
        """MicroParameters.class_fortran_names() excludes MacroParameters fields."""
        cfn = MicroParameters.class_fortran_names()
        assert "cols" not in cfn
        assert "rows" not in cfn

    def test_macro_excludes_micro_fields(self):
        """MacroParameters.class_fortran_names() excludes MicroParameters fields."""
        cfn = MacroParameters.class_fortran_names()
        assert "fiber_radius" not in cfn
        assert "micro_simulations" not in cfn

    def test_subset_of_fortran_names(self):
        """class_fortran_names() is a subset of fortran_names()."""
        all_fn = MicroParameters.fortran_names()
        cfn = MicroParameters.class_fortran_names()
        for key, val in cfn.items():
            assert all_fn[key] == val


# ===========================================================================
# Parameters.inverse_fortran_map()
# ===========================================================================


class TestInverseFortranMap:
    """inverse_fortran_map() builds {fortran_name_lower: (py, transform, cls)}."""

    def test_returns_dict(self):
        """Return type is dict."""
        result = MicroParameters.inverse_fortran_map()
        assert isinstance(result, dict)

    def test_keys_are_lowercase(self):
        """All keys in the inverse map are lowercase."""
        result = MicroParameters.inverse_fortran_map()
        for key in result:
            assert key == key.lower(), f"Key {key!r} is not lowercase"

    def test_simple_mapping(self):
        """radius maps to (fiber_radius, None, MicroParameters)."""
        result = MicroParameters.inverse_fortran_map()
        assert result["radius"] == ("fiber_radius", None, MicroParameters)

    def test_minus1_transform(self):
        """Ffree-1 maps to (empty_rows, 'minus1', MacroParameters)."""
        result = MacroParameters.inverse_fortran_map()
        assert result["ffree"] == ("empty_rows", "minus1", MacroParameters)

    def test_extra_cls_adds_entries(self):
        """extra_cls entries are included when not conflicting."""
        result = MacroParameters.inverse_fortran_map(extra_cls=MicroParameters)
        # radius is a MicroParameters field
        assert "radius" in result
        assert result["radius"][2] is MicroParameters

    def test_primary_cls_wins_conflict(self):
        """Primary cls entries override extra_cls on conflict."""
        # Both MicroParameters and MacroParameters have 'seed' Fortran name
        macro_fn = MacroParameters.class_fortran_names()
        micro_fn = MicroParameters.class_fortran_names()
        # Verify both have a 'seed' mapping (micro_seed and macro_seed)
        assert "micro_seed" in micro_fn
        assert "macro_seed" in macro_fn
        assert micro_fn["micro_seed"] == "seed|uint32"
        assert macro_fn["macro_seed"] == "seed|uint32"
        # When MacroParameters is primary, macro_seed wins
        result = MacroParameters.inverse_fortran_map(extra_cls=MicroParameters)
        assert result["seed"][0] == "macro_seed"
        assert result["seed"][2] is MacroParameters

    def test_no_extra_cls_excludes_other(self):
        """Without extra_cls, only primary class fields are included."""
        result = MicroParameters.inverse_fortran_map()
        # cols is a MacroParameters field
        macro_cfn = MacroParameters.class_fortran_names()
        def _strip_suffix(spec):
            if spec.endswith("-1"):
                return spec[:-2]
            if spec.endswith("*100"):
                return spec[:-4]
            if spec.endswith("|uint32"):
                return spec[:-7]
            return spec

        macro_fortran_lower = {
            _strip_suffix(spec).lower() for spec in macro_cfn.values()
        }
        # No macro-only Fortran names should appear
        for key in result:
            assert key not in (macro_fortran_lower - {
                k.lower() for k in MicroParameters.class_fortran_names().values()
            }) or result[key][2] is MicroParameters

    def test_tuple_has_three_elements(self):
        """Each value is a 3-tuple of (str, transform, type)."""
        result = MicroParameters.inverse_fortran_map()
        for key, val in result.items():
            assert len(val) == 3, f"Key {key!r} has {len(val)} elements"
            assert isinstance(val[0], str)
            assert val[1] in (None, "minus1", "times100", "uint32")
            assert val[2] is MicroParameters


# ===========================================================================
# Parameters.apply_fortran_transform()
# ===========================================================================


class TestApplyFortranTransform:
    """apply_fortran_transform() applies the correct numeric transform."""

    def test_none_identity(self):
        """None transform returns the value unchanged."""
        assert Parameters.apply_fortran_transform(42.0, None) == 42.0

    def test_minus1(self):
        """'minus1' subtracts 1."""
        assert Parameters.apply_fortran_transform(10.0, "minus1") == 9.0

    def test_times100(self):
        """'times100' divides by 100."""
        assert Parameters.apply_fortran_transform(50.0, "times100") == pytest.approx(0.5)

    def test_minus1_with_zero(self):
        """'minus1' on 0 returns -1."""
        assert Parameters.apply_fortran_transform(0.0, "minus1") == -1.0

    def test_times100_preserves_precision(self):
        """'times100' preserves precision for fractional values."""
        assert Parameters.apply_fortran_transform(1.0, "times100") == pytest.approx(0.01)


# ===========================================================================
# Parameters.to_quantity_or_number()
# ===========================================================================


class TestToQuantityOrNumber:
    """to_quantity_or_number() wraps values in Quantity when units exist."""

    def test_with_units_returns_quantity(self):
        """Value with known units returns a Quantity."""
        units_dict = {"fiber_radius": "microns"}
        result = Parameters.to_quantity_or_number(36.35, "fiber_radius", units_dict)
        assert isinstance(result, Quantity)
        assert result.magnitude == pytest.approx(36.35)

    def test_without_units_returns_bare_number(self):
        """Value without known units returns a bare number."""
        units_dict = {"fiber_radius": "microns"}
        result = Parameters.to_quantity_or_number(7, "nodes_in_micro_row", units_dict)
        assert result == 7
        assert not isinstance(result, Quantity)

    def test_empty_units_dict(self):
        """Empty units dict always returns bare number."""
        result = Parameters.to_quantity_or_number(100, "anything", {})
        assert result == 100
        assert not isinstance(result, Quantity)

    def test_quantity_has_correct_unit(self):
        """Returned Quantity has the expected unit."""
        units_dict = {"total_time": "sec"}
        result = Parameters.to_quantity_or_number(3600.0, "total_time", units_dict)
        assert str(result.units) == "second"
