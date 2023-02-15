# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# parameters.py
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

"""Code for holding, storing, and reading information about an experiment

This module gives a uniform way to handle the data and parameters of a given 
experiment. It contains classes to house these and make them accessible to the 
rest of the code. It also handles the storing and reading of parameters and 
data to/from disk.

Typical usage example:
    >>> # Create a new experiment
    >>> exp = Experiment('path/to/data')
    >>> param = {'override_parameter': 2.54, 'another_new_parameter': 32}
    >>> exp.initialize_params(param)
    >>> exp.to_file()
    >>> # Load an existing experiment
    >>> exp = Experiment('path/to/data', '2022_12_27_1100')
    >>> exp.read_file()
    >>> # Access a parameter
    >>> exp.params['pore_size']
    >>> # Access data
    >>> exp.data.lysis_time[4][18]
"""

import inspect
import json
import logging
import os
import pkgutil
import re
from dataclasses import asdict, dataclass, field
from typing import Any, AnyStr, List, Mapping, Tuple, Union

from ..util import DataStore, DataStatus, default_filenames, dict_to_formatted_str


class OLDExperiment:
    def __str__(self) -> str:
        """Gives a human-readable, formatted string of the current experimental
        parameters."""
        # Convert internal storage to a dictionary
        values = self.to_dict()
        # Call the formatter and return
        return dict_to_formatted_str(values)

    def to_file(self) -> None:
        """Stores the experiment parameters to disk.

        Creates or overwrites the params.json file in the experiment's data
        folder. This file will contain the current experiment parameters
        (including any micro- and macroscale parameters) in JSON format.
        """
        with open(self.os_param_file, "w") as file:
            # Convert the internal parameters to a dictionary and then use the
            # JSON module to save to disk.
            json.dump(self.to_dict(), file)


@dataclass(frozen=True)
class Parameters:
    """Contains parameters for a clot lysis scenario.

    Parameters can be accessed as attributes.
    Independent parameters should only be set at initialization.
    Dependent parameters should never be set manually, but are automatically
    calculated by internal code.

    Should only be used inside a Scenario object.
    """

    #####################################
    # Physical Parameters
    #####################################

    binding_rate: float = 0.1
    """The tPA binding rate.
    
    :Units: (micromolar*sec)^-1
    :Fortran: kon"""

    pore_size: float = 1.0135e-4
    """Pore size (distance between fibers/nodes)
    
    :Units: centimeters
    :Fortran: delx"""

    diffusion_coeff: float = 5.0e-7
    """Diffusion coefficient
    
    :Units: cm^2/s.
    :Fortran: Diff"""

    binding_sites: int = 427
    """Concentration of binding sites.
     
    :Units: micromolar
    :Fortran: bs"""

    # TODO(bpaynter): This value should derive from Microscale data
    forced_unbind: float = 0.0852
    """Fraction of times tPA was forced to unbind in microscale model.
    
    :Units: None
    :Fortran: frac_forced"""

    # TODO(bpaynter): This value should derive from Microscale data
    average_bind_time: float = 27.8
    """this is the average time a tPA molecule stays bound to fibrin. 
    For now I'm using 27.8 to be 1/0.036, the value in the absence of PLG.
    
    :Units: seconds
    :Fortran: avgwait = 1/koff"""

    # TODO(bpaynter): This value should derive from pore_size and fiber_diameter
    grid_node_distance: float = 1.0862
    """Distance from the start of one fiber to the next because distance 
    between nodes is 1.0135 micron and diameter of one fiber is 0.0727 micron.
    
    :Units: microns
    :Fortran: dist"""

    #####################################
    # Model Parameters
    #####################################

    cols: int = 93
    """The number of lattice nodes in each (horizontal) row
    
    :Units: nodes
    :Fortran: N"""

    rows: int = 121
    """The number of lattice nodes in each (vertical) column
    
    :Units: nodes
    :Fortran: F"""

    full_row: int = field(init=False)
    """Edges in a full row of nodes
    
    :Units: edges
    :Fortran: None"""

    xz_row: int = field(init=False)
    """Number of all x- and z-edges in a row
    
    :Units: edges
    :Fortran: None"""

    total_edges: int = field(init=False)
    """The total number of edges in the model
    
    :Units: edges
    :Fortran: num"""

    total_fibers: int = field(init=False)
    """The total number of fibers in the model

    :Units: fibers
    :Fortran: None"""

    empty_rows: int = 29 - 1
    """1st node in vertical direction containing fibers.
    So if first_fiber_row = 10, then rows 0-9 have no fibers, there's one more 
    row of fiber-free planar vertical edges, and then the row with index 
    'first_fiber_row' (e.g. 11th) is a full row of fibers
    
    :Units: nodes
    :Fortran: Ffree-1"""

    last_empty_edge: int = field(init=False)
    """The 1-D index of the last edge without fibrin
    
    This is probably unnecessary when using a 2-D data structure, but is kept 
    for historical reasons.
    
    :Units: edges
    :Fortran: enoFB-1"""

    total_molecules: int = 43074
    """The total number of tPA molecules:
    
        * 43074 is Colin's [tPA]=0.6 nM
        * 86148 is Colin's [tPA]=1.2 nM
        
    :Units: molecules
    :Fortran: M"""

    moving_probability: float = 0.2
    """The probability of moving.
    
    Make sure it is small enough that we've converged.
    
    :Units: None
    :Fortran: q"""

    #####################################
    # Experimental Parameters
    #####################################

    # TODO(bpaynter): This value should derive from MicroParameters
    microscale_runs: int = 50000
    """The number of independent trials run in the microscale model.
    
    :Units: trials
    :Fortran: nummicro"""

    total_trials: int = 10
    """The number of independent trials to be run
    
    :Units: trials
    :Fortran: stats"""

    total_time: int = 20 * 60
    """Total running time for model.
     
    :Units: seconds
    :Fortran: tf"""

    time_step: float = field(init=False)
    """The length of one timestep.
    
    :Units: seconds
    :Fortran: tstep"""

    total_time_steps: int = field(init=False)
    """The total number of timesteps.
    
    :Units: timesteps
    :Fortran: num_t"""

    seed: int = -2137354075
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    #####################################
    # Data Parameters
    #####################################

    input_data: List[str] = field(init=False)
    """The data (from the Microscale model) required to run the Macroscale 
    model."""

    output_data: List[str] = field(init=False)
    """The data output by the Macroscale model."""

    save_interval: int = 10
    """How often to record data from the model.
    
    :Units: sec
    :Fortran: None"""

    number_of_saves: int = field(init=False)
    """The number of times data will be saved from the model.
    
    :Units: None
    :Fortran: nplt"""

    #####################################
    # Code Parameters
    #####################################

    macro_version: str = "diffuse_into_and_along"
    """A string identifying which version of the Macroscale model is being run.
    This string was included in data filenames stored by the Fortran code."""

    log_lvl: int = logging.WARNING
    """How much debugging information to write out to the console

    :Units: None
    :Fortran: None"""

    duplicate_fortran: bool = False
    """Whether the Python code should follow the Fortran code step-by-step.
    Theoretically, with this set to "True", both sets of code will produce the 
    exact same output.
    The 
    This will impact performance negatively.

    :Units: None
    :Fortran: None"""

    processing_library: str = "numpy"
    """Which library the macroscale model should use for processing. 
    Options include
    
    * 'numpy'
    * 'cupy'
    
    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        """This method calculates the dependent parameters once the
        MacroParameters object is created. It is automatically called by the
        DataClass.__init__()"""
        # These names must be elements of the Experiment's DataStore
        object.__setattr__(
            self,
            "input_data",
            [
                "unbinding_time",  # Fortran: tsec1
                # 'leaving_time',           # Fortran: CDFtPA
                "lysis_time",  # Fortran: lysismat
                "total_lyses",  # Fortran: lenlysismat
            ],
        )
        # These names must be elements of the Experiment's DataStore
        object.__setattr__(
            self,
            "output_data",
            [
                "degradation_state",  # Fortran: degnext
                "molecule_location",
                "molecule_state",
                "save_time",  # Fortran: tsave
            ],
        )
        # A full row of the fiber grid contains a 'right', 'up', and 'out' edge
        # for each node, except the last node which contains no 'right' edge.
        object.__setattr__(self, "full_row", 3 * self.cols - 1)

        # A full row of 'right' and 'out' edges is two per node, except the
        # last node which has no 'right' edge.
        object.__setattr__(self, "xz_row", 2 * self.cols - 1)

        # The total number of edges in the grid is a full_row for each,
        # except the last row which has no 'up' edges.
        object.__setattr__(
            self, "total_edges", self.full_row * (self.rows - 1) + self.xz_row
        )

        # The total number of fibers in the grid is a full_row for each,
        # except the empty rows which have no fibers, and the last row which
        # has no 'up' edges.
        object.__setattr__(
            self,
            "total_fibers",
            self.full_row * (self.rows - self.empty_rows - 1) + self.xz_row,
        )

        # The 1-D index of the last edge in the fibrin-free region is the total
        # number of edges in the fibrin-free region -1.
        # The total rows in the fibrin-free region is equal to the (0-based)
        # index of the first fiber row.
        # The total edges in this region is one full row of edges for each row
        object.__setattr__(self, "last_empty_edge", self.full_row * self.empty_rows - 1)

        # Equation (2.4) page 25 from Bannish, et. al. 2014
        # https://doi.org/10.1093/imammb/dqs029
        object.__setattr__(
            self,
            "time_step",
            (
                self.moving_probability
                * self.pore_size**2
                / (12 * self.diffusion_coeff)
            ),
        )

        # Total timesteps is total time divided by length of one timestep
        object.__setattr__(
            self, "total_time_steps", int(self.total_time / self.time_step)
        )

        # Total saves is one for the start of each 'save_interval' plus one at
        # the end of the run.
        object.__setattr__(
            self, "number_of_saves", int(self.total_time / self.save_interval) + 1
        )

    @staticmethod
    def fortran_names():
        text = pkgutil.get_data(__name__, "parameters.py")
        pattern = re.compile(
            r"[\r\n]+^\s{4}([a-z_]+):[^\"]*\"\"\"[^\"]*:Fortran:\s([\w_]+(-1)?)(\s=[^\"]*)?\"\"\"",
            re.M,
        )
        names = {}
        matches = re.findall(pattern, text.decode("utf-8"))
        for match in matches:
            if match[1] != "None":
                names[match[0]] = match[1]
        return names

    def to_dict(self) -> Mapping[AnyStr, Any]:
        """Returns the contents of this structure as a dictionary"""
        return asdict(self)

    def __str__(self) -> str:
        """Returns a human-readable, JSON-like string of all parameters."""
        # Convert the internal parameters into one dictionary
        values = asdict(self)
        # Format the dictionary and return
        return dict_to_formatted_str(values)

    @staticmethod
    def print_default_values() -> str:
        """Returns the default parameters."""
        # Create a new MacroParameters object with the default values
        default_params = Parameters()
        # Convert to a dict, then to a formatted string, and return
        return str(default_params)
