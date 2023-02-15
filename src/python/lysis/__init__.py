#!/usr/bin/env python3
# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Brittany Bannish & Bradley Paynter
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

__version__ = "0.1"
__author__ = "Brittany Bannish, Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish & Bradley Paynter"
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]

from .edge_grid import *
from .molecule import *
from .np_macroscale import *

try:
    import cupy
    from .cp_macroscale import *
except ImportError:
    pass
