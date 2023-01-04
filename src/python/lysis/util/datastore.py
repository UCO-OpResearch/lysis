import os
from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Union

import numpy as np

from .util import dict_to_formatted_str

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

# The following are hacks to deal with loading data from Fortran into Python.
# These must be manually coded for each file
# TODO(bpaynter): Eliminate these once the whole system is converted to Python.
_data_load_args = {
                        # The total_lyses matrix needs to be converted to integer
                        'lenlysisvect.dat': {'dtype': int, 'converters': float}
                   }

_data_post_load = {
                        # The lysis_time matrix is currently stored with a full 500 rows.
                        # At the end of finite data, the rest of the column is filled with the value 6000
                        # to indicate null or infinite.
                        # Thus, we trim the lysis_time matrix to eliminate most of the rows that contain just 6000
                        # For now, I'm hardcoding this to keep 110 rows.
                        # The actual amount of data that needs to be kept is total_lyses.max()
                        #
                        # Finally we transpose the matrix for convenient access in Python (e.g. matrix[i] = row)
                        'lysismat.dat': lambda a: np.delete(a, np.s_[110:], 0).T
                  }


@unique
class Status(Flag):
    NONE = 0
    INITIALIZED = auto()
    LOADED = auto()
    SAVED = auto()


class DataStore:
    _internal_names: List[str] = ['internal_names', 'filenames', 'path', 'data', 'status']

    def __init__(self,
                 path: Union[str, bytes, os.PathLike],
                 filenames: dict[AnyStr, AnyStr] = None
                 ):
        object.__setattr__(self, '_status', {})
        object.__setattr__(self, '_path', path)
        if filenames is not None:
            object.__setattr__(self, '_filenames', filenames)
        else:
            object.__setattr__(self, '_filenames', {})
        object.__setattr__(self, '_data', {})

    def to_dict(self) -> dict[AnyStr, AnyStr]:
        return self._filenames

    def __str__(self):
        return dict_to_formatted_str(self._filenames)

    def __setattr__(self, key: AnyStr, value: Any):
        if key[0] == '_':
            raise RuntimeError(f"Don't touch my private members! ({key})")
        elif Status.INITIALIZED not in self.status(key):
            raise RuntimeError(f'{key} not initialized. Use .new() first')
        elif Status.SAVED in self.status(key):
            raise RuntimeError(f'{key} already exists (as {self._filenames[key]}). '
                               f'This module should NOT be used for modifying existing data on disk.')
        elif Status.LOADED in self.status(key):
            raise RuntimeError(f'{key} already has data in memory. Use .delete(), .overwrite(), or .append() instead.')
        else:
            self._data[key] = value
            self._set_status(key, Status.LOADED)

    def __getattr__(self, key: AnyStr) -> np.ndarray:
        if key in self._data:
            return self._data[key]
        elif Status.INITIALIZED not in self.status(key):
            raise RuntimeError(f'{key} not initialized.')
        else:
            self.load_from_disk(key)
            return self._data[key]

    def status(self, key: AnyStr) -> Status:
        if key in self._status:
            return self._status[key]
        if key not in self._filenames:
            return Status.NONE
        else:
            status = Status.INITIALIZED
            if key in self._data:
                status = status | Status.LOADED
            if os.path.isfile(os.path.join(self._path, self._filenames[key])):
                status = status | Status.SAVED
            self._status[key] = status
            return status

    def _set_status(self, key: AnyStr, status: Status):
        current_status = self.status(key)
        if current_status == Status.NONE:
            raise RuntimeError(f"This shouldn't happen. Don't try to set the status without at least a filename set.")
        self._status[key] = current_status | status

    def _unset_status(self, key: AnyStr, status: Status):
        current_status = self.status(key)
        if current_status == Status.NONE:
            raise RuntimeError(f"This shouldn't happen. Don't try to set the status without at least a filename set.")
        self._status[key] = current_status & ~status

    def new(self, key: AnyStr, filename: AnyStr):
        if key[0] == '_':
            raise RuntimeError(f"Don't touch my private members! ('{key}')")
        if key in self._internal_names:
            raise RuntimeError(f"'{key}' overrides the name of an internal attribute.")
        if Status.INITIALIZED in self.status(key):
            raise RuntimeError(f"'{key}' is already initialized with file '{self._filenames[key]}'.")
        else:
            self._filenames[key] = filename
            self._set_status(key, Status.INITIALIZED)

    def delete(self, key: AnyStr):
        if Status.SAVED in self.status(key):
            raise RuntimeError(f"'{key}' already saved to disk as '{self._filenames[key]}'. '"
                               f"This module should NOT be used for modifying existing data on disk.")
        elif key in self._data:
            self._data.pop(key)
            self._unset_status(key, Status.LOADED)
        else:
            pass

    def overwrite(self, key: AnyStr, value: Any):
        if Status.SAVED in self.status(key):
            raise RuntimeError(f"'{key}' already saved to disk as '{self._filenames[key]}'. "
                               f"This module should NOT be used for modifying existing data on disk.")
        elif key in self._data:
            self._data[key] = value
        else:
            self.__setattr__(key, value)

    def append(self, key: AnyStr, value: Any, axis: int | None = None):
        if Status.INITIALIZED not in self.status(key):
            raise RuntimeError(f'{key} not initialized.')
        elif Status.SAVED in self.status(key):
            raise RuntimeError(f'{key} already saved to disk as {self._filenames[key]}. '
                               f'This module should NOT be used for modifying existing data on disk.')
        elif Status.LOADED in self.status(key):
            self._data[key] = np.append(self._data[key], value, axis)
        else:
            self.__setattr__(key, value)

    def load_from_disk(self, key: AnyStr):
        if Status.INITIALIZED not in self.status(key):
            raise RuntimeError(f'{key} not initialized.')
        elif Status.SAVED not in self.status(key):
            raise RuntimeError(f'No file for {key} found on disk. '
                               f'({os.path.join(self._path, self._filenames[key])})')
        else:
            filename = self._filenames[key]
            args = _data_load_args.get(filename, {})
            post_load = _data_post_load.get(filename, lambda x: x)
            # If the file is stored in text format (with the extension '.dat')
            # This is the format used by the original Fortran/Matlab code
            if os.path.splitext(filename)[1] == '.dat':
                data = np.loadtxt(os.path.join(self._path, self._filenames[key]), **args)
                self._data[key] = post_load(data)
                self._set_status(key, Status.LOADED)
            # If the file is stored in NumPy format (with the extension .npy)
            # This will be the format used by the Python model once complete
            # TODO(bpaynter): Implement, consistent with the rest of the model
            elif os.path.splitext(filename)[1] == '.npy':
                raise NotImplementedError('Support for NumPy files is yet to be implemented.')
            # Other file types/extensions are not supported (or planned to be) at this time
            else:
                raise AttributeError('Non-supported file type.')

    def save_to_disk(self, key: AnyStr):
        # TODO(bpaynter): It needs to be checked whether the text files output by this method can be read by
        #   the existing Fortran and Matlab code.
        if Status.INITIALIZED not in self.status(key):
            raise RuntimeError(f'{key} not initialized.')
        elif Status.LOADED not in self.status(key):
            raise RuntimeError(f'No data for {key}.')
        elif Status.SAVED in self.status(key):
            raise RuntimeError(f'{key} already saved to disk as {self._filenames[key]}. '
                               f'This module should NOT be used for modifying existing data on disk.')
        else:
            filename = self._filenames[key]
            if os.path.splitext(filename)[1] == '.dat':
                # Save as a text file with the extension .dat
                np.savetxt(os.path.join(self._path, filename),
                           self._data[key],
                           newline=os.linesep)
            else:
                # Save as a NumPy file with the extension .npy
                np.save(os.path.join(self._path, filename))
