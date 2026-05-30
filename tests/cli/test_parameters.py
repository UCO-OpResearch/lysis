"""Tests for ``lysis parameters`` CLI command."""

import warnings
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.cli.parameters import _DEFAULT_PARAMS, _load_run_params


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


def _make_params(value="1.234", add_names=()):
    """Return a mock dict as returned by _load_run_params."""
    result = {attr_name: value for attr_name, *_ in _DEFAULT_PARAMS}
    result.update({name: value for name in add_names})
    return result


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_params():
    return _make_params()


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestParametersHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["parameters", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_markdown(self, runner):
        result = runner.invoke(cli, ["parameters", "--help"])
        assert "--markdown" in result.output

    def test_help_mentions_add_drop(self, runner):
        result = runner.invoke(cli, ["parameters", "--help"])
        assert "--add" in result.output
        assert "--drop" in result.output


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestParametersMarkdownSingleFile:
    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout_has_sections(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "### Macroscale Parameters" in result.output
        assert "### Microscale Parameters" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout_contains_param(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "pore_size" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_file(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_B.h5"
        h5.touch()
        out_file = tmp_path / "params.md"

        result = runner.invoke(
            cli, ["parameters", str(h5), "--markdown", str(out_file)]
        )

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_file_content(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_C.h5"
        h5.touch()
        out_file = tmp_path / "params.md"

        runner.invoke(cli, ["parameters", str(h5), "--markdown", str(out_file)])

        content = out_file.read_text()
        assert "## run_C" in content
        assert "### Macroscale Parameters" in content
        assert "### Microscale Parameters" in content

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_file_confirmation_message(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "params.md"

        result = runner.invoke(
            cli, ["parameters", str(h5), "--markdown", str(out_file)]
        )

        assert "Markdown written to" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        h5 = tmp_path / "run_E.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "## run_E" not in result.output


# ---------------------------------------------------------------------------
# CLI integration — directory mode
# ---------------------------------------------------------------------------


class TestParametersMarkdownDirectory:
    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        (tmp_path / "run_1.h5").touch()
        (tmp_path / "run_2.h5").touch()

        result = runner.invoke(
            cli, ["parameters", str(tmp_path), "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "### Macroscale Parameters" in result.output
        assert "### Microscale Parameters" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout_run_codes_as_columns(
        self, mock_load, runner, mock_params, tmp_path
    ):
        mock_load.return_value = (mock_params, True)
        (tmp_path / "run_X.h5").touch()
        (tmp_path / "run_Y.h5").touch()

        result = runner.invoke(
            cli, ["parameters", str(tmp_path), "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "run_X" in result.output
        assert "run_Y" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_file(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = (mock_params, True)
        (tmp_path / "run_1.h5").touch()
        out_file = tmp_path / "params.md"

        result = runner.invoke(
            cli, ["parameters", str(tmp_path), "--markdown", str(out_file)]
        )

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_with_drop(self, mock_load, runner, tmp_path):
        # Return params without the dropped key
        params = _make_params()
        del params["pore_size"]
        mock_load.return_value = (params, True)
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(
            cli,
            ["parameters", str(tmp_path), "--drop", "pore_size", "--markdown", "-"],
        )

        assert result.exit_code == 0
        assert "pore_size" not in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_with_add(self, mock_load, runner, tmp_path):
        params = _make_params(add_names=["protofibril_radius"])
        params["protofibril_radius"] = "0.0024"
        mock_load.return_value = (params, True)
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(
            cli,
            [
                "parameters",
                str(tmp_path),
                "--add",
                "protofibril_radius",
                "--markdown",
                "-",
            ],
        )

        assert result.exit_code == 0
        assert "protofibril_radius" in result.output
        assert "### Additional Parameters" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_no_markdown_flag_gives_normal_output(
        self, mock_load, runner, mock_params, tmp_path
    ):
        mock_load.return_value = (mock_params, True)
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(
            cli, ["parameters", str(tmp_path), "--no-progress"]
        )

        assert result.exit_code == 0
        assert "### Macroscale Parameters" not in result.output


# ---------------------------------------------------------------------------
# Microscale-only files (issue #22)
# ---------------------------------------------------------------------------


def _make_micro_only_run(tmp_path, run_code="micro_only", micro_simulations=42):
    """Create a microscale-only HDF5 file (no macroscale collection)."""
    from lysis.config.parameters import MicroParameters
    from lysis.dataio.datastore import DataStore

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        micro = MicroParameters(micro_simulations=micro_simulations)
    ds = DataStore.create(run_code, str(tmp_path), micro)
    ds.close()
    return run_code


class TestParametersMicroOnly:
    """A microscale-only file must not error just because no macro params exist."""

    def test_load_run_params_returns_dict(self, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        values, has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        assert values is not None

    def test_has_macro_false_for_micro_only(self, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        _values, has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        assert has_macro is False

    def test_micro_param_populated(self, tmp_path):
        run_code = _make_micro_only_run(tmp_path, micro_simulations=42)
        from rich.console import Console

        values, has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        assert values["micro_simulations"] == "42"

    def test_macro_param_is_na(self, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        values, has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        # Macro and computed params resolve to N/A, not an exception
        assert values["pore_size"] == "N/A"
        assert values["grid_node_distance"] == "N/A"

    def test_cli_single_file_exits_0(self, runner, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        h5 = tmp_path / f"{run_code}.h5"

        result = runner.invoke(cli, ["parameters", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "micro_simulations" in result.output

    def test_cli_directory_exits_0(self, runner, tmp_path):
        _make_micro_only_run(tmp_path)

        result = runner.invoke(cli, ["parameters", str(tmp_path), "--no-progress"])

        assert result.exit_code == 0

    def test_cli_drops_macro_rows(self, runner, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        h5 = tmp_path / f"{run_code}.h5"

        result = runner.invoke(
            cli, ["parameters", str(h5), "--markdown", "-"]
        )

        assert result.exit_code == 0
        # No macro values exist, so the macro section/rows are dropped entirely.
        assert "### Macroscale Parameters" not in result.output
        assert "pore_size" not in result.output
        assert "grid_node_distance" not in result.output
        # Microscale rows still render.
        assert "### Microscale Parameters" in result.output
        assert "micro_simulations" in result.output

    def test_cli_no_na_rows_remain(self, runner, tmp_path):
        run_code = _make_micro_only_run(tmp_path)
        h5 = tmp_path / f"{run_code}.h5"

        result = runner.invoke(cli, ["parameters", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "N/A" not in result.output

    def test_new_micro_params_real_values(self, tmp_path):
        """The five micro Fortran-CLI inputs added in #91 show real values."""
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        values, _has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        assert values["nodes_in_micro_row"] == "7"
        assert values["snap_proportion"] == "0.6667"
        assert values["unbind_rate_PLi"] == "57.60"
        assert values["activation_rate_PLG"] == "0.100"
        assert values["exposure_rate_binding_site"] == "5.00"

    def test_micro_params_broad_coverage_not_na(self, tmp_path):
        """Every microscale entry in the default table resolves (issue #91)."""
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        values, _has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )

        micro_attrs = [attr for attr, src, *_ in _DEFAULT_PARAMS if src == "micro"]
        na = [a for a in micro_attrs if values[a] == "N/A"]
        assert na == [], f"micro params unexpectedly N/A: {na}"

    def test_add_micro_fallback_not_in_defaults(self, tmp_path):
        """``--add`` of a micro attr absent from the defaults still resolves.

        ``protofibril_radius`` is a microscale parameter that is not part of
        ``_DEFAULT_PARAMS``; pulling it via ``--add`` on a micro-only file
        exercises the ``_get_raw`` microscale fallback that issue #91 targets.
        """
        run_code = _make_micro_only_run(tmp_path)
        from rich.console import Console

        values, _has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, ["protofibril_radius"], Console()
        )

        assert values["protofibril_radius"] != "N/A"


# ---------------------------------------------------------------------------
# Macro-spec dropping (canonical source)
# ---------------------------------------------------------------------------


class TestDropMacroSpecs:
    """``_drop_macro_specs`` keys off the canonical MacroParameters fields."""

    def test_macro_param_names_from_dataclass(self):
        from lysis.cli.parameters import _macro_param_names

        names = _macro_param_names()
        # Representative macroscale-only fields
        assert "pore_size" in names
        assert "rows" in names
        assert "macro_simulations" in names
        # Microscale fields and the nested container are excluded
        assert "fiber_radius" not in names
        assert "micro_simulations" not in names
        assert "micro_params" not in names

    def test_drops_macro_and_computed_keeps_micro(self):
        from lysis.cli.parameters import _drop_macro_specs

        kept = {p[0] for p in _drop_macro_specs(_DEFAULT_PARAMS)}
        # Macro and computed rows are gone
        assert "pore_size" not in kept
        assert "grid_node_distance" not in kept  # computed
        # Every microscale row survives, so a missing micro param still shows N/A
        assert "fiber_radius" in kept
        assert "micro_simulations" in kept
        assert "micro_seed" in kept


# ---------------------------------------------------------------------------
# Non-default value highlighting (rich terminal only)
# ---------------------------------------------------------------------------


def _render_rich(df):
    """Render a parameters DataFrame to an ANSI string via a forced terminal."""
    import io

    from rich.console import Console

    from lysis.tools.display import params_df_to_rich

    buf = io.StringIO()
    console = Console(
        file=buf, force_terminal=True, width=120, color_system="standard"
    )
    console.print(params_df_to_rich(df))
    return buf.getvalue()


class TestNonDefaultHighlight:
    """Values differing from the model defaults are highlighted in rich output."""

    def _flags_for(self, tmp_path, **micro_kwargs):
        from rich.console import Console

        from lysis.cli.parameters import _default_formatted, _nondefault_flags

        run_code = _make_micro_only_run(tmp_path, **micro_kwargs)
        values, _has_macro = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )
        flags = _nondefault_flags(values, _default_formatted(_DEFAULT_PARAMS, []))
        return run_code, values, flags

    def _df_for(self, tmp_path, **micro_kwargs):
        from lysis.analysis.summary import parameters_table
        from lysis.cli.parameters import _drop_macro_specs
        from lysis.config.parameters import MacroParameters

        run_code, values, flags = self._flags_for(tmp_path, **micro_kwargs)
        specs = _drop_macro_specs(_DEFAULT_PARAMS)
        return parameters_table(
            {run_code: values},
            specs,
            [],
            MacroParameters.units(),
            [run_code],
            nondefault_by_run={run_code: flags},
        )

    def test_nondefault_flags_mark_changed_only(self, tmp_path):
        # micro_simulations default is 50000; 42 is non-default.
        _run_code, _values, flags = self._flags_for(tmp_path, micro_simulations=42)

        assert flags["micro_simulations"] is True
        # Left at their defaults:
        assert flags["bind_rate_tPA"] is False
        assert flags["snap_proportion"] is False
        assert flags["nodes_in_micro_row"] is False

    def test_na_cells_never_flagged(self, tmp_path):
        """Macro params (absent here, rendered N/A) are never marked non-default."""
        _run_code, values, flags = self._flags_for(tmp_path)

        assert values["pore_size"] == "N/A"
        assert flags["pore_size"] is False

    def test_table_attaches_aligned_nondefault_frame(self, tmp_path):
        df = self._df_for(tmp_path, micro_simulations=42)

        flags_df = df.attrs.get("nondefault")
        assert flags_df is not None
        assert list(flags_df.index) == list(df.index)
        assert list(flags_df.columns) == list(df.columns)

    def test_rich_render_colors_nondefault_cell(self, tmp_path):
        out = _render_rich(self._df_for(tmp_path, micro_simulations=42))

        # The non-default micro_simulations row carries the yellow escape; a
        # left-at-default row (bind_rate_tPA) does not.  Assert on \x1b[33m
        # specifically: every row's bold "Parameter" cell already emits \x1b[,
        # so a bare \x1b[ check would pass even with highlighting removed.
        nondefault_line = next(
            l for l in out.splitlines() if "micro_simulations" in l
        )
        default_line = next(l for l in out.splitlines() if "bind_rate_tPA" in l)
        assert "\x1b[33m" in nondefault_line
        assert "\x1b[33m" not in default_line

    def test_no_attrs_means_no_color(self, tmp_path):
        """Without a nondefault frame, the rich table is rendered plain."""
        from rich.console import Console

        from lysis.analysis.summary import parameters_table
        from lysis.cli.parameters import _drop_macro_specs
        from lysis.config.parameters import MacroParameters

        run_code = _make_micro_only_run(tmp_path, micro_simulations=42)
        values, _ = _load_run_params(
            str(tmp_path), run_code, _DEFAULT_PARAMS, [], Console()
        )
        specs = _drop_macro_specs(_DEFAULT_PARAMS)
        df = parameters_table(
            {run_code: values}, specs, [], MacroParameters.units(), [run_code]
        )  # no nondefault_by_run

        assert "\x1b[33m" not in _render_rich(df)

    def test_markdown_output_is_unstyled(self, runner, tmp_path):
        run_code = _make_micro_only_run(tmp_path, micro_simulations=42)
        h5 = tmp_path / f"{run_code}.h5"

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "\x1b[" not in result.output
