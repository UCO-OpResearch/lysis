# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# data_file.py
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

"""Structures and methods for dealing with individual data files

"""

import os

from dataclasses import dataclass, replace
from enum import Enum, unique
from shutil import copy2
from typing import Any, AnyStr, Callable, Dict, Self

import numpy as np
import pandas as pd


@unique
class DataType(Enum):
    FORTRAN = 0
    MATLAB = 1
    NUMPY = 2
    PANDAS = 3


load_command = {
    DataType.FORTRAN: np.fromfile,
    DataType.MATLAB: np.loadtxt,
    DataType.NUMPY: np.load,
    DataType.PANDAS: pd.read_json,
}

extension = {
    DataType.FORTRAN: ".dat",
    DataType.MATLAB: ".dat",
    DataType.NUMPY: ".npy",
    DataType.PANDAS: ".json",
}

save_command = {
    DataType.NUMPY: np.save,
    DataType.PANDAS: lambda x, y: y.to_json(x),
}


@dataclass
class DataFile:
    filename: AnyStr = None
    path: AnyStr = None
    data_type: DataType = None
    post_load: Callable = None
    load_args: Dict = None
    contents: Any = None

    def can_save(self) -> bool:
        valid = True
        valid = valid and self.filename is not None
        valid = valid and self.path is not None
        valid = valid and self.contents is not None
        valid = valid and self.data_type in DataType
        valid = (
            valid and os.path.splitext(self.filename)[1] == extension[self.data_type]
        )
        valid = valid and os.path.isdir(self.path)
        valid = valid and not os.path.isfile(os.path.join(self.path, self.filename))
        return valid

    def can_load(self) -> bool:
        valid = True
        valid = valid and self.filename is not None
        valid = valid and self.path is not None
        valid = valid and self.data_type in DataType
        valid = (
            valid and os.path.splitext(self.filename)[1] == extension[self.data_type]
        )
        valid = valid and self.contents is None
        valid = valid and os.path.isfile(os.path.join(self.path, self.filename))
        return valid

    def load(self):
        if not self.can_load():
            raise RuntimeError(f"Cannot load file {self.filename}.")
        if self.load_args is None:
            load_args = {}
        else:
            load_args = self.load_args
        self.contents = load_command[self.data_type](
            os.path.join(self.path, self.filename), **load_args
        )
        if self.post_load is not None:
            self.contents = self.post_load(self.contents)

    def save(self):
        if not self.can_save():
            raise RuntimeError(f"Cannot save file {self.filename}.")
        if self.data_type not in save_command:
            raise RuntimeError(f"Cannot save data of type {self.data_type}.")
        save_command[self.data_type](
            os.path.join(self.path, self.filename), self.contents
        )

    def as_type(self, data_type: DataType) -> Self:
        if self.data_type == data_type:
            return self
        out = DataFile(data_type=data_type)
        out.filename = os.path.splitext(self.filename)[0] + extension[data_type]
        out.path = self.path
        out.contents = self.contents
        out.post_load = self.post_load
        out.load_args = self.load_args
        return out

    def move(self, new_path):
        if os.path.isfile(os.path.join(self.path, self.filename)):
            os.rename(
                os.path.join(self.path, self.filename),
                os.path.join(new_path, self.filename),
            )
        self.path = new_path

    def copy(self, new_path) -> Self:
        if os.path.isfile(os.path.join(self.path, self.filename)):
            copy2(os.path.join(self.path, self.filename), new_path)
        return replace(self, path=new_path)
