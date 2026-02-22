"""Unit tests for lysis.config.paramcheck.

Tests cover the six public functions:
    load_micro_params, load_macro_params      -- strict parameter loading
    parse_micro_log, parse_macro_log          -- Fortran log parsing
    verify_micro_params, verify_macro_params  -- log-vs-params verification

All log-parsing tests use inline text constants (no real data files).
"""

import textwrap
import warnings

import pytest

from lysis.config.constants import Q_
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.config.paramcheck import (
    load_macro_params,
    load_micro_params,
    parse_macro_log,
    parse_micro_log,
    verify_macro_params,
    verify_micro_params,
)


# ---------------------------------------------------------------------------
# Inline Fortran log content (no dependency on actual data files)
# ---------------------------------------------------------------------------

# Representative micro log through the first stats= line.
# All values match MicroParameters defaults except:
#   nodes=13  (nodes_in_micro_row=7 by default)
#   seed=2133256963  (micro_seed=0 by default)
# Includes "Setting" lines from command-line argument parsing at the top,
# followed by the standard key=value parameter echo.
MICRO_LOG_CONTENT = textwrap.dedent("""\
 number of command arguments =            8
 command arg            1 = nodes
 command arg            2 =           13
 Setting nodes =           13
 command arg            3 = simulations
 command arg            4 =        50000
 Setting simulations =        50000
 command arg            5 = seed
 command arg            6 =   2133256963
 Setting seed =   2133256963
 command arg            7 = outFileCode
 command arg            8 = _PLG2_tPA01_TB-xiii
 Setting outFileCode = _PLG2_tPA01_TB-xiii
 command line processed
  filetype=binary
  seed=  2133256963
 nodes=          13
 KdtPAnoplg=  0.360000000000000
 KdtPAyesplg=  2.000000000000000E-002
 KdPLGnicked=   2.20000000000000
 KdPLGintact=   38.0000000000000
 simulations=       50000
 data/2024-04-16-1920/lysis__PLG2_tPA01_TB-xiii.dat
  kncat=   5.00000000000000
  kapcat=  0.100000000000000
  ktPAon=  0.100000000000000
  kaoff10=  3.600000000000000E-002
  kaoff12=  2.000000000000000E-003
  kplioff=   57.6000000000000
  kplgoffnick=  0.220000000000000
  kplgon=  0.100000000000000
  kplgoff=   3.80000000000000
  freeplg=   2.00000000000000
  kdeg=   5.00000000000000
  stats=        1000
  stats=        2000
""")

# Representative macro log params section.
# Values: N=19, F=184, Ffree=113 (→ empty_rows=112), M=21105, seed=-725030,
# q=0.2, delx=5.34e-4 cm (pore_size), frac_forced=0.074101.
# Integer dependent params included: num=10285 (total_edges), enoFB=6272 (empty_edges).
# Excluded: avgwait (floating-point mismatch), tstep (floating-point mismatch),
#           kon (requires overrides), bs (unit mismatch with stored value).
MACRO_LOG_CONTENT = textwrap.dedent("""\
 command line processed
  N=          19
  F=         184
  Ffree=         113
  num=       10285
  M=                 21105
  seed=     -725030
 enoFB=        6272
  run number=                     1
  q=  0.200000000000000
  delx=  5.340000000000000E-004
  tf=           0
  num_t=  0.000000000000000E+000
  frac_forced=  7.410100000000000E-002
After     10. sec, 0 fibers are degraded
""")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_log(content, tmp_path, filename="test_log.txt"):
    """Write inline content to a temp file and return its Path."""
    path = tmp_path / filename
    path.write_text(content)
    return path


def _make_micro():
    """Return default MicroParameters (no warnings context needed)."""
    return MicroParameters()


def _make_macro(micro=None, **kwargs):
    """Return MacroParameters, suppressing any unit-assumption warnings."""
    if micro is None:
        micro = _make_micro()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return MacroParameters(micro_params=micro, **kwargs)


# ---------------------------------------------------------------------------
# TestLoadMicroParams
# ---------------------------------------------------------------------------


class TestLoadMicroParams:
    """Tests for load_micro_params(): strict loading with dependency checks."""

    def test_full_valid_data(self, capsys):
        """Round-trip to_basedict() → load_micro_params() succeeds with correct values."""
        micro = _make_micro()
        bd = micro.to_basedict()
        loaded = load_micro_params(bd)
        capsys.readouterr()  # discard debug print from parse_from_basedict

        assert loaded.micro_simulations == micro.micro_simulations
        assert loaded.micro_seed == micro.micro_seed
        assert loaded.nodes_in_micro_row == micro.nodes_in_micro_row

    def test_missing_independent_raises(self, capsys):
        """Deleting one independent param raises ValueError naming it."""
        micro = _make_micro()
        bd = micro.to_basedict()
        del bd["fiber_radius"]

        with pytest.raises(ValueError, match="fiber_radius"):
            load_micro_params(bd)
        capsys.readouterr()

    def test_multiple_missing_raises(self, capsys):
        """Deleting two independent params raises ValueError naming both."""
        micro = _make_micro()
        bd = micro.to_basedict()
        del bd["fiber_radius"]
        del bd["nodes_in_micro_row"]

        with pytest.raises(ValueError, match="fiber_radius") as exc:
            load_micro_params(bd)
        assert "nodes_in_micro_row" in str(exc.value)
        capsys.readouterr()

    def test_override_fills_missing(self, capsys):
        """A direct-value override supplies a missing independent param."""
        micro = _make_micro()
        bd = micro.to_basedict()
        del bd["fiber_radius"]

        loaded = load_micro_params(bd, overrides={"fiber_radius": "36.35 nanometer"})
        capsys.readouterr()

        assert loaded.fiber_radius.to("nm").magnitude == pytest.approx(36.35)

    def test_dependent_mismatch_raises(self, capsys):
        """A tampered stored dependent param raises ValueError naming the field."""
        micro = _make_micro()
        bd = micro.to_basedict()
        # protofibril_radius = 2 × fibrinogen_radius = 0.0024 µm; tamper it
        bd["protofibril_radius"] = "0.001 micron"

        with pytest.raises(ValueError, match="protofibril_radius"):
            load_micro_params(bd)
        capsys.readouterr()

    def test_tolerance_accepted(self, capsys):
        """A stored dependent value within relative tolerance does not raise."""
        micro = _make_micro()
        bd = micro.to_basedict()
        # Introduce a sub-tolerance relative error (1e-10 << default 1e-9)
        exact = micro.protofibril_radius.to("microns").magnitude  # 0.0024
        bd["protofibril_radius"] = f"{exact * (1 + 1e-10)} micron"

        load_micro_params(bd)  # must not raise
        capsys.readouterr()


# ---------------------------------------------------------------------------
# TestLoadMacroParams
# ---------------------------------------------------------------------------


class TestLoadMacroParams:
    """Tests for load_macro_params(): strict loading for MacroParameters."""

    def test_full_valid_data(self, capsys):
        """Round-trip to_basedict() → load_macro_params() preserves key fields."""
        micro = _make_micro()
        macro = _make_macro(micro)
        bd = macro.to_basedict()
        loaded = load_macro_params(bd, micro)
        capsys.readouterr()

        assert loaded.cols == macro.cols
        assert loaded.rows == macro.rows
        assert loaded.empty_rows == macro.empty_rows

    def test_missing_independent_raises(self, capsys):
        """Deleting one independent macro param raises ValueError naming it."""
        micro = _make_micro()
        macro = _make_macro(micro)
        bd = macro.to_basedict()
        del bd["cols"]

        with pytest.raises(ValueError, match="cols"):
            load_macro_params(bd, micro)
        capsys.readouterr()

    def test_override_fills_missing(self, capsys):
        """A direct-value override supplies a missing macro independent param."""
        micro = _make_micro()
        macro = _make_macro(micro)
        bd = macro.to_basedict()
        del bd["cols"]

        loaded = load_macro_params(bd, micro, overrides={"cols": 93})
        capsys.readouterr()

        assert loaded.cols == 93

    def test_dependent_mismatch_raises(self, capsys):
        """A tampered stored dependent param raises ValueError naming the field."""
        micro = _make_micro()
        macro = _make_macro(micro)
        bd = macro.to_basedict()
        # full_row = 3*cols - 1 = 278; tamper it
        bd["full_row"] = 999

        with pytest.raises(ValueError, match="full_row"):
            load_macro_params(bd, micro)
        capsys.readouterr()

    def test_requires_micro_params_not_none(self):
        """Passing micro_params=None raises ValueError immediately."""
        micro = _make_micro()
        macro = _make_macro(micro)
        bd = macro.to_basedict()

        with pytest.raises(ValueError, match="micro_params"):
            load_macro_params(bd, None)


# ---------------------------------------------------------------------------
# TestParseMicroLog
# ---------------------------------------------------------------------------


class TestParseMicroLog:
    """Tests for parse_micro_log(): Fortran micro log → Python param dict."""

    def test_parses_dissociation_constants(self, tmp_path):
        """KdtPAnoplg= maps to diss_const_tPA_woPLG with correct value and units."""
        path = _write_log(MICRO_LOG_CONTENT, tmp_path)
        result = parse_micro_log(path)

        assert "diss_const_tPA_woPLG" in result
        parsed = result["diss_const_tPA_woPLG"]
        assert parsed.to("micromolar").magnitude == pytest.approx(0.36)

    def test_parses_kplgon(self, tmp_path):
        """kplgon= maps to bind_rate_PLG (case-insensitive key matching)."""
        path = _write_log(MICRO_LOG_CONTENT, tmp_path)
        result = parse_micro_log(path)

        assert "bind_rate_PLG" in result
        assert result["bind_rate_PLG"].magnitude == pytest.approx(0.1)

    def test_stops_at_stats_line(self, tmp_path):
        """Lines at and after the first stats= line are not parsed."""
        content = textwrap.dedent("""\
             KdtPAnoplg=  0.360000000000000
              stats=        1000
              KdtPAyesplg=  0.999
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "diss_const_tPA_woPLG" in result
        # diss_const_tPA_wPLG appears AFTER stats=; must not be parsed
        assert "diss_const_tPA_wPLG" not in result

    def test_unknown_name_raises(self, tmp_path):
        """A numeric-valued Fortran name absent from fortran_names() raises ValueError."""
        content = " runs=       50000\n stats=1\n"
        path = _write_log(content, tmp_path)

        with pytest.raises(ValueError, match="runs"):
            parse_micro_log(path)

    def test_setting_lines_parsed(self, tmp_path):
        """'Setting key = value' lines from command-line parsing are extracted."""
        content = textwrap.dedent("""\
             Setting nodes =           13
             Setting simulations =        50000
             Setting seed =   2133256963
             Setting outFileCode = _PLG2_tPA01_TB-xiii
              stats=        1000
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "nodes_in_micro_row" in result
        assert result["nodes_in_micro_row"] == pytest.approx(13)
        assert "micro_simulations" in result
        assert result["micro_simulations"] == pytest.approx(50000)
        assert "micro_seed" in result
        assert result["micro_seed"] == pytest.approx(2133256963)
        # String-valued Setting lines (outFileCode) are silently skipped
        assert "outFileCode" not in result

    def test_setting_overwritten_by_later_kv(self, tmp_path):
        """A later key=value line overwrites an earlier Setting line."""
        content = textwrap.dedent("""\
             Setting nodes =           99
             nodes=          13
              stats=        1000
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        # The later nodes=13 should overwrite Setting nodes=99
        assert result["nodes_in_micro_row"] == pytest.approx(13)

    def test_alias_override_resolves_unknown(self, tmp_path):
        """overrides={"micro_simulations": "runs"} maps runs= to micro_simulations."""
        content = " runs=       50000\n stats=1\n"
        path = _write_log(content, tmp_path)

        result = parse_micro_log(path, overrides={"micro_simulations": "runs"})

        assert "micro_simulations" in result
        assert result["micro_simulations"] == pytest.approx(50000)


# ---------------------------------------------------------------------------
# TestParseMacroLog
# ---------------------------------------------------------------------------


class TestParseMacroLog:
    """Tests for parse_macro_log(): Fortran macro log → Python param dict."""

    def test_parses_grid_params(self, tmp_path):
        """N=, F=, M= map to cols, rows, total_molecules with correct values."""
        path = _write_log(MACRO_LOG_CONTENT, tmp_path)
        result = parse_macro_log(path)

        assert result["cols"] == pytest.approx(19)
        assert result["rows"] == pytest.approx(184)
        assert result["total_molecules"] == pytest.approx(21105)

    def test_ffree_transform(self, tmp_path):
        """Ffree=113 → empty_rows=112  (Python = Fortran − 1)."""
        path = _write_log(MACRO_LOG_CONTENT, tmp_path)
        result = parse_macro_log(path)

        assert "empty_rows" in result
        assert result["empty_rows"] == pytest.approx(112)

    def test_stops_at_after_line(self, tmp_path):
        """Lines beginning with 'After ' are not parsed."""
        content = textwrap.dedent("""\
             N=          19
             After     10. sec, 0 fibers are degraded
             F=         999
        """)
        path = _write_log(content, tmp_path)
        result = parse_macro_log(path)

        assert result["cols"] == pytest.approx(19)
        # rows (F=) appears after the After line; must not be parsed
        assert "rows" not in result

    def test_cross_class_params_resolved(self, tmp_path):
        """bs= in macro log resolves to MicroParameters.binding_sites."""
        content = " bs=         42\nAfter done\n"
        path = _write_log(content, tmp_path)
        result = parse_macro_log(path)

        assert "binding_sites" in result
        assert result["binding_sites"].magnitude == pytest.approx(42)


# ---------------------------------------------------------------------------
# TestVerifyMicroParams
# ---------------------------------------------------------------------------


class TestVerifyMicroParams:
    """Tests for verify_micro_params(): log values matched against params object."""

    def _micro_matching_log(self):
        """MicroParameters consistent with MICRO_LOG_CONTENT."""
        return MicroParameters(
            nodes_in_micro_row=13,
            micro_simulations=50000,
            micro_seed=2133256963,
        )

    def test_matching_params_no_error(self, tmp_path):
        """Params consistent with log raise no exception."""
        micro = self._micro_matching_log()
        path = _write_log(MICRO_LOG_CONTENT, tmp_path)
        verify_micro_params(micro, path)  # must not raise

    def test_mismatch_raises(self, tmp_path):
        """A log value that differs from stored params raises ValueError."""
        # Log has KdtPAnoplg=0.36; params have diss_const_tPA_woPLG=0.36 µM but
        # we swap to a clearly wrong value in a fresh log.
        content = " KdtPAnoplg=  0.500000000000000\n stats=1\n"
        path = _write_log(content, tmp_path)
        micro = _make_micro()  # diss_const_tPA_woPLG = 0.36 µM (mismatch)

        with pytest.raises(ValueError, match="diss_const_tPA_woPLG"):
            verify_micro_params(micro, path)

    def test_tolerance_accepted(self, tmp_path):
        """A log value within relative tolerance does not raise."""
        micro = _make_micro()
        # Introduce a relative error well below 1e-9
        exact = micro.diss_const_tPA_woPLG.to("micromolar").magnitude  # 0.36
        perturbed = exact * (1 + 1e-11)
        content = f" KdtPAnoplg=  {perturbed:.15E}\n stats=1\n"
        path = _write_log(content, tmp_path)

        verify_micro_params(micro, path)  # must not raise

    def test_direct_value_override(self, tmp_path):
        """A direct-value override replaces the stored attribute for comparison."""
        # Log has non-default KdtPAnoplg
        content = " KdtPAnoplg=  0.500000000000000\n stats=1\n"
        path = _write_log(content, tmp_path)
        micro = _make_micro()  # default 0.36 µM → mismatch without override

        # Without override: should raise (0.5 != 0.36)
        with pytest.raises(ValueError):
            verify_micro_params(micro, path)

        # With override: 0.5 µM is used as the expected value → matches log
        verify_micro_params(
            micro, path, overrides={"diss_const_tPA_woPLG": "0.5 micromolar"}
        )


# ---------------------------------------------------------------------------
# TestVerifyMacroParams
# ---------------------------------------------------------------------------


class TestVerifyMacroParams:
    """Tests for verify_macro_params(): macro log values matched against params."""

    def _macro_matching_log(self):
        """MacroParameters consistent with MACRO_LOG_CONTENT."""
        micro = MicroParameters(nodes_in_micro_row=13, micro_seed=2133256963)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return MacroParameters(
                micro_params=micro,
                cols=19,
                rows=184,
                empty_rows=112,          # Ffree=113 → empty_rows=112
                total_molecules=21105,
                macro_seed=-725030,
                moving_probability=0.2,
                pore_size=Q_("5.34e-4 cm"),
                forced_unbind=0.074101,
                total_time=Q_("0 sec"),  # tf=0 in log
            )

    def test_matching_params_no_error(self, tmp_path):
        """MacroParameters consistent with MACRO_LOG_CONTENT raise no exception."""
        macro = self._macro_matching_log()
        path = _write_log(MACRO_LOG_CONTENT, tmp_path)
        verify_macro_params(macro, path)  # must not raise

    def test_mismatch_raises(self, tmp_path):
        """A log value differing from stored params raises ValueError."""
        # Log has delx=5.34e-4; params use wrong pore_size
        content = " delx=  5.340000000000000E-004\nAfter done\n"
        path = _write_log(content, tmp_path)
        wrong_macro = _make_macro(pore_size=Q_("1.0 cm"))  # clearly wrong

        with pytest.raises(ValueError, match="pore_size"):
            verify_macro_params(wrong_macro, path)

    def test_ffree_verified_correctly(self, tmp_path):
        """Ffree=N in log is verified against stored empty_rows = N−1."""
        # Ffree=5 → empty_rows should be 4
        content = " Ffree=          5\nAfter done\n"
        path = _write_log(content, tmp_path)

        # Correct: empty_rows=4 matches Ffree=5
        correct = _make_macro(cols=10, rows=10, empty_rows=4)
        verify_macro_params(correct, path)  # must not raise

        # Wrong: empty_rows=5 does NOT match Ffree=5 (that would require empty_rows=4)
        wrong = _make_macro(cols=10, rows=10, empty_rows=5)
        with pytest.raises(ValueError, match="empty_rows"):
            verify_macro_params(wrong, path)
