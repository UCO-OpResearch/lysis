"""Shared fixtures for lysis.dataio unit tests."""

import os
import pathlib

import pytest
import numpy as np

from lysis.config.constants import CONST
from lysis.dataio.dataspec import DataSetSpec, DataCollectionSpec


# ----------------------------------------------------------------------
# Suppress the src/lysis/ dirty-tree gate during the test session.
#
# The provenance feature added in ``provenance-stamps`` errors out when
# ``src/lysis/`` has uncommitted changes — but developers iterating in
# the source tree expect every test to pass regardless of working-tree
# state.  Setting the env-var override globally means the tests behave
# as if the tree were clean; the gate itself is exercised directly in
# ``tests/cli/test_provenance_gates.py`` where the env var is cleared.
# ----------------------------------------------------------------------
@pytest.fixture(autouse=True, scope="session")
def _allow_provenance_gates_during_tests():
    os.environ.setdefault(CONST.LYSIS_ALLOW_DIRTY_ENV, "1")
    os.environ.setdefault(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV, "1")
    yield


@pytest.fixture
def sample_params():
    """Small BaseParamsType dict for tests that need parameter resolution."""
    return {
        "micro_params": {
            "micro_simulations": 100,
        },
        "macro_params": {
            "rows": 5,
            "cols": 3,
            "total_edges": 57,
            "total_molecules": 10,
            "empty_rows": 0,
            "full_row": 8,
            "pore_size": 1.0,
        },
    }


@pytest.fixture
def text_spec():
    """DataSetSpec for a simple text file dataset."""
    return DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
        dtype=np.float64,
        data_location="data{file_code}.dat",
        shape=(-1,),
    )


@pytest.fixture
def binary_spec():
    """DataSetSpec for a simple binary file dataset."""
    return DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
        dtype=np.int32,
        data_location="data{file_code}.dat",
        shape=(3, 4),
    )


@pytest.fixture
def json_spec():
    """DataSetSpec for a JSON parameter file."""
    return DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_JSON,
        dtype=dict,
        data_location="params{file_code}.json",
        shape=(-1,),
    )


@pytest.fixture
def hdf5_dataset_spec():
    """DataSetSpec for an HDF5 dataset."""
    return DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
        dtype=np.float64,
        data_location="data/arr",
        shape=(-1,),
    )


@pytest.fixture
def hdf5_attr_spec():
    """DataSetSpec for HDF5 attributes."""
    return DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
        dtype=dict,
        data_location="test_data",
        shape=(-1,),
    )


@pytest.fixture
def combined_collection_spec(json_spec, text_spec):
    """DataCollectionSpec with simulations_combined=True."""
    return DataCollectionSpec(
        simulations_combined=True,
        params=json_spec,
        data={"arr": text_spec},
    )


@pytest.fixture
def per_sim_collection_spec(json_spec):
    """DataCollectionSpec with simulations_combined=False and {sim:02} paths."""
    sim_spec = DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
        dtype=np.float64,
        data_location="{sim:02}/data{file_code}.dat",
        shape=(-1,),
    )
    return DataCollectionSpec(
        simulations_combined=False,
        params=json_spec,
        data={"arr": sim_spec},
    )


# ──────────────────────────────────────────────────────────────────────
# Real data fixtures
# ──────────────────────────────────────────────────────────────────────

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_FIXTURE_DIR = _REPO_ROOT / "tests" / "fixtures" / "fortran_sample"
_FULL_DATA_DIR = _REPO_ROOT / "data" / "2026-02-18-1723"


@pytest.fixture(scope="session")
def fortran_sample_path():
    """Path to the committed truncated Fortran data fixture.

    Returns the path to tests/fixtures/fortran_sample/, which contains
    a truncated subset of real Fortran simulation output.
    """
    if not _FIXTURE_DIR.exists():
        pytest.skip("Truncated fixture not found at tests/fixtures/fortran_sample/")
    return str(_FIXTURE_DIR)


@pytest.fixture(scope="session")
def full_data_path():
    """Path to the full (~665MB) Fortran data directory (local only).

    Tests using this fixture should be marked with @pytest.mark.real_data
    so they are skipped in CI.
    """
    if not _FULL_DATA_DIR.exists():
        pytest.skip("Full dataset not found at data/2026-02-18-1723/")
    return str(_FULL_DATA_DIR)


_FIXTURE_V190_DIR = _REPO_ROOT / "tests" / "fixtures" / "fortran_v190_sample"
_FULL_DATA_V190_DIR = _REPO_ROOT / "data" / "2026-02-28-1907"


@pytest.fixture(scope="session")
def fortran_v190_sample_path():
    """Path to the committed truncated v1.90.0 Fortran data fixture.

    Returns the path to tests/fixtures/fortran_v190_sample/, which contains
    a truncated subset of real Fortran v1.90.0 simulation output (3 snapshots,
    simulation 00 only).
    """
    if not _FIXTURE_V190_DIR.exists():
        pytest.skip("Truncated v1.90.0 fixture not found at tests/fixtures/fortran_v190_sample/")
    return str(_FIXTURE_V190_DIR)


@pytest.fixture(scope="session")
def full_data_v190_path():
    """Path to the full (~190MB) v1.90.0 Fortran data directory (local only).

    Tests using this fixture should be marked with @pytest.mark.real_data
    so they are skipped in CI.
    """
    if not _FULL_DATA_V190_DIR.exists():
        pytest.skip("Full v1.90.0 dataset not found at data/2026-02-28-1907/")
    return str(_FULL_DATA_V190_DIR)
