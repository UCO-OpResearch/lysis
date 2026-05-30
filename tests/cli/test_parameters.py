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
