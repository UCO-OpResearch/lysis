"""Shared fixtures for lysis.data unit tests."""

import pytest
import numpy as np

from lysis.config.constants import CONST
from lysis.data.dataspec import DataSetSpec, DataCollectionSpec


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
