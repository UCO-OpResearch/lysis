"""Read legacy seed attributes without rewriting their HDF5 files."""

import h5py
import numpy as np

from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.dataio.datastore import DataStore


def test_run_loads_legacy_float_seed(tmp_path):
    with DataStore.create("legacy-run", str(tmp_path), MicroParameters()):
        pass
    path = tmp_path / "legacy-run.h5"
    with h5py.File(path, "a") as data:
        data.attrs["converted_from"] = "v1.95.0"
        data["micro_data"].attrs["micro_seed"] = np.float64(1093845299.0)
    stored_bytes = path.read_bytes()

    run = Run(tmp_path, "legacy-run")
    run.load_params_from_hdf5()

    assert run.micro_params.seed_as_int() == 1093845299
    assert run.micro_params.micro_seed == "1093845299"
    assert path.read_bytes() == stored_bytes
