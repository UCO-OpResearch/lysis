# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# fortran_simulation.py
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

"""Simulation object wrapper for Fortran code."""

import inspect
import os
import subprocess

from dataclasses import asdict, dataclass
from typing import AnyStr

from experiment import Model, Scenario, Simulation


class FortranSimulation(Simulation):
    def __init__(self, model: Model, scenario: Scenario, index: int = 0):
        super().__init__(model, scenario, index)

    def run(self):
        pass
