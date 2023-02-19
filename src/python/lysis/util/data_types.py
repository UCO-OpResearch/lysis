# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# data_types.py
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

"""Implementations of the `DataType` class

"""
from typing import AnyStr, Any

import numpy as np

from .data_file import DataType


pandas = True
try:
    import pandas as pd
except ImportError:
    pandas = False
    pd = None

cupy = True
try:
    import cupy as cp
except ImportError:
    cupy = False
    cp = None


__all__ = ["DataType", "NumPyData", "FortranData", "MatlabData"]
if pandas:
    __all__.append("PandasJsonData")
if cupy:
    __all__.append("CuPyData")


class NumPyData(DataType):
    @staticmethod
    def save(filepath: AnyStr, data: np.ndarray, *args, **kwargs):
        return np.save(filepath, data, *args, **kwargs)

    @staticmethod
    def load(filepath: AnyStr, *args, **kwargs) -> np.ndarray:
        return np.load(filepath, *args, **kwargs)

    @property
    def extension(self) -> AnyStr:
        return ".npy"

    @property
    def name(self) -> AnyStr:
        return "NumPy"


class FortranData(DataType):
    """Reading Fortran binary files as NumPy Arrays.

    Writing not implemented."""

    @staticmethod
    def save(filepath: AnyStr, data: Any, *args, **kwargs):
        """This method is not implemented and will raise a
        `NotImplementedError` if called."""
        raise NotImplementedError("Saving as Fortran Binary is not supported.")

    @staticmethod
    def load(filepath: AnyStr, *args, **kwargs):
        return np.fromfile(filepath, *args, **kwargs)

    @property
    def extension(self) -> AnyStr:
        return ".dat"

    @property
    def name(self) -> AnyStr:
        return "Fortran Binary"


class MatlabData(DataType):
    """Reading Matlab data files as NumPy Arrays.

    Writing not implemented."""

    @staticmethod
    def save(filepath: AnyStr, data: Any, *args, **kwargs):
        raise NotImplementedError("Saving as Matlab is not supported.")

    @staticmethod
    def load(filepath: AnyStr, *args, **kwargs):
        return np.loadtxt(filepath, *args, **kwargs)

    @property
    def extension(self) -> AnyStr:
        return ".dat"

    @property
    def name(self) -> AnyStr:
        return "Matlab"


class PandasJsonData(DataType):
    @staticmethod
    def save(filepath: AnyStr, data: pd.DataFrame, *args, **kwargs):
        return data.to_json(filepath, *args, **kwargs)

    @staticmethod
    def load(filepath: AnyStr, *args, **kwargs):
        return pd.read_json(filepath, *args, **kwargs)

    @property
    def extension(self) -> AnyStr:
        return ".json"

    @property
    def name(self) -> AnyStr:
        return "Pandas JSON"


class CuPyData(DataType):
    @staticmethod
    def save(filepath: AnyStr, data: Any, *args, **kwargs):
        return np.save(filepath, cp.asnumpy(data), *args, **kwargs)

    @staticmethod
    def load(filepath: AnyStr, *args, **kwargs):
        return cp.array(np.load(filepath, *args, **kwargs))

    @property
    def extension(self) -> AnyStr:
        return ".npy"

    @property
    def name(self) -> AnyStr:
        return "CuPy"
