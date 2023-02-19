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
from shutil import copy2
from typing import Any, AnyStr, Callable, Concatenate, Dict, ParamSpec, Self

P = ParamSpec("P")


@dataclass(frozen=True)
class DataType:
    name: AnyStr
    extension: AnyStr
    save: Callable[Concatenate[AnyStr, Any, P], Any] = None
    load: Callable[Concatenate[AnyStr, P], Any] = None


@dataclass
class DataFile:
    filename: AnyStr = None
    load_args: Dict = None
    post_load: Callable = None
    save_args: Dict = None
    pre_save: Callable = None
    data_type: DataType = None
    contents: Any = None

    def can_save(self) -> bool:
        valid = True
        valid = valid and self.filename is not None
        valid = valid and self.contents is not None
        valid = valid and self.data_type is not None
        valid = valid and os.path.splitext(self.filename)[1] == self.data_type.extension
        valid = valid and os.path.isdir(os.path.dirname(self.filename))
        valid = valid and not os.path.isfile(self.filename)
        return valid

    def can_load(self) -> bool:
        valid = True
        valid = valid and self.filename is not None
        valid = valid and self.data_type is not None
        valid = valid and os.path.splitext(self.filename)[1] == self.data_type.extension
        valid = valid and self.contents is None
        valid = valid and os.path.isfile(self.filename)
        return valid

    def load(self):
        if not self.can_load():
            raise RuntimeError(f"Cannot load file {self.filename}.")
        if self.load_args is None:
            load_args = {}
        else:
            load_args = self.load_args
        self.contents = self.data_type.load(self.filename, **load_args)
        if self.post_load is not None:
            self.contents = self.post_load(self.contents)

    def save(self):
        if not self.can_save():
            raise RuntimeError(f"Cannot save file {self.filename}.")
        if self.data_type.save is None:
            raise NotImplementedError(
                f"Saving is not implemented for data type {self.data_type.name}."
            )
        if self.save_args is None:
            save_args = {}
        else:
            save_args = self.save_args
        out = self.contents
        if self.pre_save is not None:
            out = self.pre_save(out)
        self.data_type.save(self.filename, out, **save_args)

    def as_type(self, data_type: DataType) -> Self:
        if self.data_type == data_type:
            return self
        out = DataFile(data_type=data_type)
        out.filename = os.path.splitext(self.filename)[0] + data_type.extension
        out.contents = self.contents
        return out

    def move(self, new_path):
        if not os.path.isdir(new_path):
            raise RuntimeError(f"{new_path} is not a valid location.")
        new_filename = os.path.join(new_path, os.path.basename(self.filename))
        if os.path.isfile(self.filename):
            os.rename(
                self.filename,
                new_filename,
            )
        self.filename = new_filename

    def copy(self, new_path) -> Self:
        if not os.path.isdir(new_path):
            raise RuntimeError(f"{new_path} is not a valid location.")
        new_filename = os.path.join(new_path, os.path.basename(self.filename))
        if os.path.isfile(self.filename):
            copy2(self.filename, new_path)
        return replace(self, filename=new_filename)
