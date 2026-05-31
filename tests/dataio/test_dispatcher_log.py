"""Tests for dispatcher-log capture (issue #55), Stage 1: reader + spec + convert.

Covers:

* ``_read_file_text`` now passes ``comments=None`` so ``#``-containing log lines
  are preserved (the worker ``.out`` may contain ``#``; numeric data files do not).
* The ``micro_dispatcher_log`` / ``macro_dispatcher_log`` DataSetSpecs exist at
  every spec version with the expected storage type and ``optional=True``.
* Identity converters are registered for the dispatcher logs in both the explicit
  v1.99.0<->v2.0.0 tables and (via spec inheritance) the v1.95.0<->v1.99.0
  auto-comprehension tables, so cross-version conversion never KeyErrors.
"""

import numpy as np

from lysis.config.constants import CONST
from lysis.dataio.dataspec import DataSetSpec, dataspec
from lysis.dataio.dataconvert import data_converters
from lysis.dataio.fileops import _read_file_text

# Pull the null delimiter straight from the loaded spec so the test never has to
# embed a raw NUL byte in source.
_DISP_DELIM = (
    dataspec["v1.99.0"]["microscale_out"].data["micro_dispatcher_log"].delimiter
)


class TestReadFileTextKeepsHashLines:
    """``comments=None`` keeps ``#`` lines that the old default would drop."""

    def _disp_spec(self):
        return DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=str,
            data_location="dlog{file_code}.out",
            delimiter=_DISP_DELIM,
        )

    def test_hash_lines_preserved(self, tmp_path):
        (tmp_path / "dlog.out").write_text(
            "first line\n# a full comment line\nmid # inline part\nlast line\n"
        )
        rows = [str(x) for x in _read_file_text(str(tmp_path), self._disp_spec())]
        # A line that begins with '#' survives (old comments='#' would drop it).
        assert any(r.lstrip().startswith("#") for r in rows)
        # Content after an inline '#' is not truncated away.
        assert any("inline part" in r for r in rows)

    def test_set_x_trace_line_preserved(self, tmp_path):
        # set -x traces are '+'-prefixed; ensure they round-trip one row per line.
        (tmp_path / "dlog.out").write_text("+ mkdir -p /tmp/x\n+ echo done\n")
        rows = [str(x) for x in _read_file_text(str(tmp_path), self._disp_spec())]
        assert rows == ["+ mkdir -p /tmp/x", "+ echo done"]


class TestDispatcherSpecsPresent:
    """The dispatcher specs exist at every version with the right shape."""

    def test_source_specs(self):
        for coll, key in [
            ("microscale_out", "micro_dispatcher_log"),
            ("macroscale_out", "macro_dispatcher_log"),
        ]:
            for ver in ("v1.99.0", "v1.95.0", "v1.90.0"):
                spec = dataspec[ver][coll].data[key]
                assert spec.dataset_storage_type == CONST.DATASET_STORAGE_TYPE.FILE_TEXT
                assert spec.optional is True
                assert spec.data_location.endswith(".out")

    def test_target_specs(self):
        micro = dataspec["v2.0.0"]["microscale_out"].data["micro_dispatcher_log"]
        macro = dataspec["v2.0.0"]["macroscale_out"].data["macro_dispatcher_log"]
        assert micro.data_location == "log_files/dispatcher/micro_dispatcher_log"
        assert (
            macro.data_location
            == "log_files/dispatcher/macro_dispatcher_log__sim_{sim:02}"
        )
        for spec in (micro, macro):
            assert spec.dataset_storage_type == CONST.DATASET_STORAGE_TYPE.HDF5_DATASET
            assert spec.optional is True


class TestDispatcherConverters:
    """Identity converters are registered in every relevant direction."""

    def test_explicit_v199_v200_tables(self):
        for table in [("v1.99.0", "v2.0.0"), ("v2.0.0", "v1.99.0")]:
            conv = data_converters[table]
            assert "micro_dispatcher_log" in conv
            assert "macro_dispatcher_log" in conv
            payload = {"micro_dispatcher_log": np.array(["a", "b"], dtype=object)}
            assert (
                conv["micro_dispatcher_log"](payload) is payload["micro_dispatcher_log"]
            )

    def test_inherited_v195_v199_tables(self):
        # v1.95.0/v1.90.0 inherit the specs via deepcopy, so the auto-comprehension
        # converter tables include the dispatcher keys -> no KeyError mid-chain.
        for table in [("v1.95.0", "v1.99.0"), ("v1.99.0", "v1.95.0")]:
            conv = data_converters[table]
            assert "micro_dispatcher_log" in conv
            assert "macro_dispatcher_log" in conv
