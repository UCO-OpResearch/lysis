"""Tests for ``lysis parameters`` CLI command."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.cli.parameters import _DEFAULT_PARAMS


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
# Wrappers around the new public API (replaces removed private helpers)
# ---------------------------------------------------------------------------


def _params_to_markdown(param_specs, add_names, natural_units, rows, ordered_codes,
                        single_file_code=None):
    from lysis.analysis.summary import parameters_table
    from lysis.tools.display import params_df_to_markdown
    df = parameters_table(rows, param_specs, add_names, natural_units, ordered_codes)
    return params_df_to_markdown(df, single_code=single_file_code)


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
# Markdown helpers (unit tests)
# ---------------------------------------------------------------------------


class TestMdTableParameters:
    def test_header_pipe_format(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["Parameter", "Run1"], [["pore_size", "1.0135"]])
        assert out.splitlines()[0] == "| Parameter | Run1 |"

    def test_separator_row(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["A", "B", "C"], [])
        assert out.splitlines()[1] == "| --- | --- | --- |"

    def test_data_row_present(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["P", "V"], [["pore_size (microns)", "1.0135"]])
        assert "| pore_size (microns) | 1.0135 |" in out


class TestParamsToMarkdown:
    """Unit tests for the parameters markdown output (via parameters_table + params_df_to_markdown)."""

    def _natural_units(self):
        from lysis.config.parameters import MacroParameters
        return MacroParameters.units()

    def test_single_file_heading(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"], single_file_code="run_A",
        )
        assert md.startswith("## run_A")

    def test_macroscale_section_header(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"],
        )
        assert "### Macroscale Parameters" in md

    def test_microscale_section_header(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"],
        )
        assert "### Microscale Parameters" in md

    def test_macro_param_appears_in_macro_section(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"],
        )
        macro_section = md.split("### Macroscale Parameters")[1].split("### Microscale")[0]
        assert "pore_size" in macro_section

    def test_micro_param_appears_in_micro_section(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"],
        )
        micro_section = md.split("### Microscale Parameters")[1]
        assert "bind_rate_tPA" in micro_section

    def test_run_code_as_column(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_X": mock_params}, ["run_X"],
        )
        assert "run_X" in md

    def test_multiple_run_codes_as_columns(self):
        rows = {"run_X": _make_params("1.0"), "run_Y": _make_params("2.0")}
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            rows, ["run_X", "run_Y"],
        )
        header_lines = [l for l in md.splitlines() if "run_X" in l]
        assert any("run_Y" in l for l in header_lines)

    def test_add_names_create_additional_section(self, mock_params):
        extra_params = {**mock_params, "protofibril_radius": "0.0024"}
        md = _params_to_markdown(
            _DEFAULT_PARAMS, ["protofibril_radius"], self._natural_units(),
            {"run_A": extra_params}, ["run_A"],
        )
        assert "### Additional Parameters" in md
        assert "protofibril_radius" in md

    def test_no_add_names_no_additional_section(self, mock_params):
        md = _params_to_markdown(
            _DEFAULT_PARAMS, [], self._natural_units(),
            {"run_A": mock_params}, ["run_A"],
        )
        assert "### Additional Parameters" not in md

    def test_drop_removes_param(self):
        filtered = [p for p in _DEFAULT_PARAMS if p[0] != "pore_size"]
        rows = {"run_A": _make_params()}
        del rows["run_A"]["pore_size"]

        md = _params_to_markdown(
            filtered, [], self._natural_units(), rows, ["run_A"],
        )
        assert "pore_size" not in md


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestParametersMarkdownSingleFile:
    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = mock_params
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout_has_sections(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = mock_params
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "### Macroscale Parameters" in result.output
        assert "### Microscale Parameters" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_stdout_contains_param(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = mock_params
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["parameters", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "pore_size" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_markdown_to_file(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = mock_params
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
        mock_load.return_value = mock_params
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
        mock_load.return_value = mock_params
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "params.md"

        result = runner.invoke(
            cli, ["parameters", str(h5), "--markdown", str(out_file)]
        )

        assert "Markdown written to" in result.output

    @patch("lysis.cli.parameters._load_run_params")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_params, tmp_path):
        mock_load.return_value = mock_params
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
        mock_load.return_value = mock_params
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
        mock_load.return_value = mock_params
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
        mock_load.return_value = mock_params
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
        mock_load.return_value = params
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
        mock_load.return_value = params
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
        mock_load.return_value = mock_params
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(
            cli, ["parameters", str(tmp_path), "--no-progress"]
        )

        assert result.exit_code == 0
        assert "### Macroscale Parameters" not in result.output
