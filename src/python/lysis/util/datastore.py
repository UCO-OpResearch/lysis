import os

from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Mapping, Union

import numpy as np
import h5py

from .constants import CONST

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


tpa_molecule_status = h5py.enum_dtype(
    {i.name: i.value for i in CONST.MOL_STATUS}, basetype="u1"
)

tpa_bind_event_type = np.dtype(
    [
        ("Simulation Time Elapsed", np.float64),
        ("tPA Molecule Index", np.int64),
        ("Molecule New Status", tpa_molecule_status),
        ("Grid Location Index", np.int16),
    ]
)

fiber_degrade_event_type = np.dtype(
    [
        ("Simulation Time Elapsed", np.float64),
        ("Grid Location Index", np.int32),
        ("Fiber New Degrade Time", np.float64),
    ]
)


@unique
class DataStatus(Flag):
    NONE = 0
    INITIALIZED = auto()
    LOADED = auto()
    SAVED = auto()
    FILLED = auto()


def h5_tree(val: h5py.Dataset, pre: AnyStr = ""):
    """Recursively prints the tree of an HDF5 file's contents.

    Copied from https://stackoverflow.com/questions/61133916/is-there-in-python-a-single-function-that-shows-the-full-structure-of-a-hdf5-fi

    Args:
        val: The item in the HDF5 to print
        pre: The current indentation
    """
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                print(pre + "└── " + key)
                h5_tree(val, pre + "    ")
            else:
                try:
                    print(pre + "└── " + key + " (%d)" % len(val))
                except TypeError:
                    print(pre + "└── " + key + " (scalar)")
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + "├── " + key)
                h5_tree(val, pre + "│   ")
            else:
                try:
                    print(pre + "├── " + key + " (%d)" % len(val))
                except TypeError:
                    print(pre + "├── " + key + " (scalar)")


class DataStore:
    """
    A simple data store class that allows for storing and retrieving data.
    """

    _internal_names: List[str] = [
        "internal_names",
        "run_code",
        "path",
        "data",
        "tables",
        "views",
        "status",
        "mode",
    ]

    def __init__(self, run_code: AnyStr, path: AnyStr, mode: str = "r"):
        """
        Initialize the DataStore with an optional run parameter.

        :param run: The experimental run associated with the data store.
        :param mode: Mode in which to open the file.

            ``'r'``
                Read only

            ``'w'``
                Write (NOTE: This will overwrite any current file contents)

            ``'a'``
                Append. Will allow the addition of new data only.
                Will NOT allow the insertion of data if some already exists.
        """
        object.__setattr__(self, "_run_code", run_code)
        object.__setattr__(self, "_path", path)
        object.__setattr__(self, "_mode", mode)
        object.__setattr__(
            self,
            "_status",
            {
                "self": DataStatus.INITIALIZED,
                "micro": DataStatus.NONE,
                "macro": DataStatus.NONE,
            },
        )

        object.__setattr__(self, "_tables", [])
        object.__setattr__(self, "_views", [])

        # Initialize HDF5 file
        h5py.get_config().track_order = True
        object.__setattr__(
            self,
            "_data",
            h5py.File(
                os.path.join(self._path, f"{self._run_code}.h5"),
                self._mode,
            ),
        )

        if "micro_data" in self._data:
            self._status["micro"] = DataStatus.INITIALIZED
        else:
            self._status["micro"] = DataStatus.NONE

        if "macro_data" in self._data:
            self._status["macro"] = DataStatus.INITIALIZED
        else:
            self._status["macro"] = DataStatus.NONE

        # TODO: Add code here to check if there is actually data stored in this HDF5
        #       and if it matches the current data specification

    def __str__(self):
        """Print the datastore in human-readable format."""
        return h5_tree(self._data)

    def import_fortran_micro_data(
        self, filecode: AnyStr = None, data_version: AnyStr = "current"
    ):
        """
        Import data from a Fortran Microscale run into the HDF5 storage.

        :param filecode: The file code associated with the Microscale run being imported.
            This file code should include any leading underscores, but NOT the file extension.
        """
        if self._mode == "r":
            raise os.UnsupportedOperation("Data is open in read-only mode.")
        if self._mode == "a" and self._status["micro"] == DataStatus.INITIALIZED:
            raise os.UnsupportedOperation("Existing data cannot be overwritten.")
        micro_data = self._data.create_group("micro_data")

    def import_fortran_macro_data(self, filecode=None):
        """
        Import data from a Fortran Macroscale run into the HDF5 storage.

        :param filecode: The file code associated with the Macroscale run being imported.
            This file code should include any leading underscores, but NOT the file extension.
        """
        pass

    def export_fortran_micro_data(self, filecode=None):
        """
        Export Microscale data from the HDF5 storage to disk for use by a Fortran Macroscale run.

        See the "Macro to Micro files" section of the Data Specification for more information.

        :param filecode: The file code associated with the Microscale run being exported.
            This file code should include any leading underscores, but NOT the file extension.
        """
        # TODO: Add code to check if microscale data exists
        # Get the microscale data
        micro_data = self._data["micro_data"]
        # Get the number of microscale runs and set the dimensions of the bins so that we get 100 bins
        bin_size = micro_data["pli_first_time"].size // 100
        # tPAleave is the CDF of the tPA leaving time distribution.
        # This is really just a list of edgepoints from the bins for tPA leaving time
        # These bins are evenly distributed along the interval [0, 1]
        tPAleave = np.append(np.arange(0, 1, 0.01), [1.0])
        np.savetxt(os.path.join(self._path, f"tPAleave{filecode}.dat"), tPAleave)

        # The remaining data will be arranged into 100 bins
        # according to the time tPA left the simulation.
        # Get the sorted ordering of the tPA leaving times
        indices = micro_data["tpa_leaving_time"][:].argsort()
        # Find the tPA leaving times for the edges of each bin.
        tsectPA = np.append(
            [0], micro_data["tpa_leaving_time"][:][indices[bin_size - 1 :: bin_size]]
        )
        np.savetxt(os.path.join(self._path, f"tsectPA{filecode}.dat"), tsectPA)

        # Identify which simulations had the fiber fully degraded
        lysis_complete = micro_data["fiber_degraded"][:]
        # Read in the fiber degradation times for all simulations
        lysis_time = micro_data["sim_final_time"][:]
        # If full degradation did NOT occur,
        # this matrix currently contains the ending time of the simulation.
        # Replace these times with an 'infinity' marker of 6,000 seconds
        lysis_time[~lysis_complete] = 6_000
        # Rearrange the matrix so that each row contains the lysis times for simulations
        # corresponding to the matching bin in the ``tsectPA`` vector.
        # Then sort the rows (bins) individually by lysis time.
        # Finally, transpose the matrix so that the bins are arranged in columns.
        lysismat = np.stack(
            [
                np.sort(lysis_time[indices[i * bin_size : (i + 1) * bin_size]])
                for i in range(100)
            ]
        ).T
        np.savetxt(os.path.join(self._path, f"lysismat{filecode}.dat"), lysismat)

        # Find the location of the first '6000' entry in each column of the ``lysismat`` matrix
        # Then convert to 1-indexing.
        lenlysisvect = lysismat.argmax(axis=1) + 1
        np.savetxt(
            os.path.join(self._path, f"lenlysisvect{filecode}.dat"), lenlysisvect
        )

    def __getattr__(self, key: AnyStr) -> np.ndarray:
        """
        Get data from the data store.
        :param key: The table name to retrieve.
        :return: The attribute value.
        """
        # Check if the attribute exists in the HDF5 file
        # If it does, return the value
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def __setattr__(self, key: AnyStr, value: np.ndarray):
        """
        Set data in the data store.
        :param key: The table name to set.
        :param value: The value to set.
        """
        pass

    def status(self, key: AnyStr) -> DataStatus:
        """
        Get the status of the data in the data store.
        :param key: The table name to check.
        :return: The status of the data.
        """
        # Check if the table exists in the HDF5 file
        # If it does, return the status
        # If it doesn't, return DataStatus.NONE
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def _set_status(self, key: AnyStr, status: DataStatus):
        """
        Set the status of the data in the data store.
        :param key: The table name to set the status for.
        :param status: The status to set.
        """
        # Check if the table exists in the HDF5 file
        # If it does, set the status
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def _unset_status(self, key: AnyStr, status: DataStatus):
        """
        Unset the status of the data in the data store.
        :param key: The table name to unset the status for.
        :param status: The status to unset.
        """
        # Check if the attribute exists in the HDF5 file
        # If it does, unset the status
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def delete(self, key: AnyStr):
        """
        Delete data from the data store.
        :param key: The table name to delete.
        """
        # Check if the table exists in the HDF5 file
        # If it does, delete it
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def overwrite(self, key: AnyStr, value: np.ndarray):
        """
        Overwrite data in the data store.
        :param key: The table name to overwrite.
        :param value: The value to overwrite with.
        """
        # Check if the table exists in the HDF5 file
        # If it does, overwrite it
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass

    def append(self, key: AnyStr, value: np.ndarray, axis: int | None = None):
        """
        Append data to the data store.
        :param key: The table name to append to.
        :param value: The value to append.
        :param axis: The axis to append along (optional).
        """
        # Check if the table exists in the HDF5 file
        # If it does, append the array
        # If it doesn't, raise an AttributeError
        # This is a placeholder implementation
        # Replace with actual logic to access HDF5 file
        pass
