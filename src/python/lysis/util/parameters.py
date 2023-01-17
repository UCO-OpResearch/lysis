"""Code for holding, storing, and reading information about an experiment

This module gives a uniform way to handle the data and parameters of a given experiment.
It contains classes to house these and make them accessible to the rest of the code.
It also handles the storing and reading of parameters and data to/from disk.

Typical usage example:
    >>> # Create a new experiment
    >>> exp = Experiment('path/to/data')
    >>> param = {'override_parameter': 2.54, 'another_new_parameter': 32}
    >>> exp.initialize_macro_param(param)
    >>> exp.to_file()
    >>> # Load an existing experiment
    >>> exp = Experiment('path/to/data', '2022_12_27_1100')
    >>> exp.read_file()
    >>> # Access a parameter
    >>> exp.macro_params['pore_size']
    >>> # Access data
    >>> exp.data.lysis_time[4][18]
"""

import inspect
import json
import logging
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, List, Mapping, Tuple, Union

from .constants import default_filenames
from .datastore import DataStore, DataStatus
from .util import dict_to_formatted_str

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class Experiment(object):
    """Houses all information about a given experimental run.

    This object contains:

    * Data location
    * Experiment parameters
    * Experiment data

    It includes methods for

    * Initializing with default parameters
    * Reading parameters from disk
    * Saving parameters to disk
    * Reading input data from disk
    * Writing result data to disk

    Args:
        data_root: The path of the folder containing datasets
        experiment_code: The code number of the experiment.
            This will be the name of the folder containing the data specific to this experiment.
            This should be a date and time in 'YYYY-MM-DD-hhmm' format
            If no code is given, one will be generated from the current date and time.

    Attributes:
        experiment_code (str): The code number of the experiment.
        os_path (str): The path to the folder containing this experiment's data
        macro_params (DataClass): A dictionary of


    Raises:
        RuntimeError: An invalid data folder was given.
    """
    def __init__(self,
                 data_root: Union[str, bytes, os.PathLike],
                 experiment_code: str = None):
        # Check if the data folder path is valid
        if not os.path.isdir(data_root):
            raise RuntimeError('Data folder not found.')
        self.os_data_root = data_root
        # If no experiment code was given, create a new one from the current date and time.
        if experiment_code is None:
            self.experiment_code = datetime.now().strftime("%Y-%m-%d-%H%M")
        else:
            self.experiment_code = experiment_code

        # Generate the path to the experiment folder and the parameters file
        self.os_path = os.path.join(data_root, str(self.experiment_code))
        self.os_param_file = os.path.join(self.os_path, 'params.json')

        # TODO(bpaynter): Check if the parameters are already stored.
        #                 Don't allow parameters to be changed once stored.
        # Initialize the internal storage as empty
        self.micro_params: Mapping[str, Any] | None = None
        self.macro_params = None
        self.data = DataStore(self.os_path, default_filenames)

    def __str__(self) -> str:
        """Gives a human-readable, formatted string of the current experimental parameters."""
        # Convert internal storage to a dictionary
        values = self.to_dict()
        # Call the formatter and return
        return dict_to_formatted_str(values)

    def initialize_macro_param(self, params: Mapping[str, Any] = None) -> None:
        """Creates the parameters for the Macroscale model.

        Parameters are set to the default values unless new values are passed in the params dictionary.

        This method is essentially a wrapper for the MacroParameters constructor.

        Args:
            params: A dictionary of parameters that differ from the default values.

                For example,
                    >>> {'binding_rate': 10, 'pore_size': 3,}
        """
        if params is not None:
            self.macro_params = MacroParameters(**params)
        else:
            self.macro_params = MacroParameters()

        for data in self.macro_params.data_required:
            if DataStatus.INITIALIZED not in self.data.status(data):
                raise RuntimeError(f"'{data}' not initialized in DataStore.")

    def to_dict(self) -> dict:
        """Returns the internally stored data as a dictionary.

        Does not include system-specific information like paths.
        """
        # Initialize a dictionary of the appropriate parameters
        output = {'experiment_code': self.experiment_code,
                  'data_filenames': None,
                  'micro_params': None,
                  'macro_params': None
                  }
        # Get the data filenames from the DataStore
        if self.data is not None:
            output['data_filenames'] = self.data.to_dict()
        # Convert the Macroscale parameters to a dictionary
        if self.macro_params is not None:
            output["macro_params"] = asdict(self.macro_params)
        # Convert the Microscale parameters to a dictionary
        # if self.micro_params is not None:
        #     output["micro_params"] = asdict(self.macro_params)
        return output
    
    def to_file(self) -> None:
        """Stores the experiment parameters to disk.

        Creates or overwrites the params.json file in the experiment's data folder.
        This file will contain the current experiment parameters (including any micro- and macroscale parameters)
        in JSON format.
        """
        with open(self.os_param_file, 'w') as file:
            # Convert the internal parameters to a dictionary and then use the JSON module to save to disk.
            json.dump(self.to_dict(), file)

    def read_file(self) -> None:
        """Load the experiment parameters from disk.

        Raises:
            RuntimeError: No parameter file is available for this experiment.
        """
        # Determine whether the parameter file exists for this experiment
        if not os.path.isfile(self.os_param_file):
            raise RuntimeError('Experiment parameter file not found.')
        # Open the file
        with open(self.os_param_file, 'r') as file:
            # Use the JSON library to read in the parameters as a dictionary
            params = json.load(file)
            data_filenames = params.pop('data_filenames', None)
            if data_filenames is not None:
                self.data = DataStore(self.os_path, data_filenames)
            # Remove the Macroscale parameters from the dictionary (if it exists)
            # and create a new MacroParameters object using its values
            macro_params = params.pop('macro_params', None)
            if macro_params is not None:
                # We are checking here to make sure that saved, dependent parameters
                # don't get passed to the MacroParameters constructor

                # Find the parameters needed to initialize a new MacroParameters object
                sig = inspect.signature(MacroParameters)
                # Get the keys we read from the JSON
                for key in list(macro_params.keys()):
                    # If that key is not needed, then toss it
                    if key not in sig.parameters:
                        macro_params.pop(key)
                # Now unpack whatever is left in the dict and pass it to the constructor
                self.macro_params = MacroParameters(**macro_params)
            else:
                # If there were no parameters in the file, then we leave the object null.
                self.macro_params = None


@dataclass(frozen=True)
class MicroParameters:
    """This will contain the parameters for the Microscale model.
    This will need to be implemented and the Fortran microscale code modified to work with it.

    Once implemented, many Macroscale parameters will need to be redefined so that they derive from the appropriate
    Microscale parameters.
    """
    # TODO(bpaynter): Needs to be implemented. This implementation should include:
    #                   * Standard Microscale parameters
    #                   * The inclusion of standard sets of Microscale parameters (i.e., Q2, CaseA-D, etc.)
    #                   * The modification of the Microscale Fortran code to read/write JSON parameters
    #                   * The modification of the MacroParameters class to derive the appropriate parameters
    #                       from the MicroParameters class
    pass


@dataclass(frozen=True)
class MacroParameters:
    """Contains parameters for the Macroscale model.

    Parameters can be accessed as attributes.
    Independent parameters should only be set at initialization.
    Dependent parameters should never be set manually, but are automatically calculated by internal code.

    Should only be used inside an Experiment object.


    Example:
        >>> # Initialize using the default values
        >>> macro_params_default = MacroParameters()
        >>> # Get parameter value
        >>> macro_params_default.pore_size
        1.0135e-4
        >>> # Initialize overriding some default values
        >>> p = {'binding_rate': 10, 'pore_size': 3}
        >>> macro_params_override = MacroParameters(**p)
    """
    #####################################
    # Physical Parameters
    #####################################

    # TODO(bpaynter): This value should derive from MicroParameters
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

    # TODO(bpaynter): This value should derive from MicroParameters
    binding_sites: float = 4.27e+2
    """Concentration of binding sites.
     
    :Units: micromolar
    :Fortran: bs"""

    # TODO(bpaynter): This value should derive from MicroParameters
    forced_unbind: float = 0.0852
    """Fraction of times tPA was forced to unbind in microscale model.
    
    :Units: None
    :Fortran: frac_forced"""

    # TODO(bpaynter): This value should derive from MicroParameters
    average_bind_time: float = 27.8
    """this is the average time a tPA molecule stays bound to fibrin. 
    For now I'm using 27.8 to be 1/0.036, the value in the absence of PLG.
    
    :Units: seconds
    :Fortran: avgwait = 1/koff"""

    # TODO(bpaynter): This value should derive from pore_size and MicroParameters['fiber_diameter']
    grid_node_distance: float = 1.0862
    """Distance from the start of one fiber to the next 
    because distance between nodes is 1.0135 micron and diameter of 1 fiber is 0.0727 micron.
    
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
    So if first_fiber_row = 10, then rows 0-9 have no fibers, there's one more row of fiber-free planar vertical edges,
    and then the row with index 'first_fiber_row' (e.g. 11th) is a full row of fibers
    
    :Units: nodes
    :Fortran: Ffree-1"""

    last_empty_edge: int = field(init=False)
    """The 1-D index of the last edge without fibrin
    
    This is probably unnecessary when using a 2-D data structure, but is kept for historical reasons.
    
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

    seed: int = (-2137354075)
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    state: Tuple[int, int, int, int] = (129281,
                                        362436069,
                                        123456789,
                                        None)
    """State for the random number generator.
    
    :Units: None
    :Fortran: state"""

    #####################################
    # Data Parameters
    #####################################

    data_required: List[str] = field(init=False)
    """The data (from the Microscale model) required to run the Macroscale model."""

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

    macro_version: str = 'diffuse_into_and_along'
    """A string identifying which version of the Macroscale model is being run.
    This string was included in data filenames stored by the Fortran code."""

    log_lvl: int = logging.WARNING
    """How much debugging information to write out to the console

    :Units: None
    :Fortran: None"""

    duplicate_fortran: bool = False
    """Whether the Python code should follow the Fortran code step-by-step.
    Theoretically, with this set to "True", both sets of code will produce the exact same output.
    This will impact performance negatively.
    This currently does nothing.

    :Units: None
    :Fortran: None"""

    processing_library: str = 'numpy'
    """Which library the macroscale model should use for processing. 
    Options include
    
    * 'numpy'
    * 'cupy'
    
    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        """This method calculates the dependent parameters once the MacroParameters object is created.
        It is automatically called by the DataClass.__init__()"""
        # These names must be elements of the Experiment's DataStore
        object.__setattr__(self, 'data_required', [
                                                    'unbinding_time',           # Fortran: tsec1
                                                    # 'leaving_time',           # Fortran: CDFtPA
                                                    'lysis_time',               # Fortran: lysismat
                                                    'total_lyses',              # Fortran: lenlysismat
                                                   ])
        # A full row of the fiber grid contains a 'right', 'up', and 'out' edge for each node,
        # except the last node which contains no 'right' edge.
        object.__setattr__(self, 'full_row', 3 * self.cols - 1)

        # A full row of 'right' and 'out' edges is two per node, except the last node which has no 'right' edge.
        object.__setattr__(self, 'xz_row', 2 * self.cols - 1)

        # The total number of edges in the grid is a full_row for each,
        # except the last row which has no 'up' edges.
        object.__setattr__(self, 'total_edges', self.full_row * (self.rows - 1) + self.xz_row)

        # The total number of fibers in the grid is a full_row for each,
        # except the empty rows which have no fibers, and the last row which has no 'up' edges.
        object.__setattr__(self, 'total_fibers', self.full_row * (self.rows - self.empty_rows - 1) + self.xz_row)

        # The 1-D index of the last edge in the fibrin-free region is the total number of edges in the
        # fibrin-free region -1.
        # The total rows in the fibrin-free region is equal to the (0-based) index of the first fiber row
        # The total edges in this region is one full row of edges for each row
        object.__setattr__(self, 'last_empty_edge', self.full_row * self.empty_rows - 1)

        # Equation (2.4) page 25 from Bannish, et. al. 2014
        # https://doi.org/10.1093/imammb/dqs029
        object.__setattr__(self, 'time_step', (self.moving_probability
                                               * self.pore_size ** 2
                                               / (12 * self.diffusion_coeff)))

        # Total timesteps is total time divided by length of one timestep
        object.__setattr__(self, 'total_time_steps', int(self.total_time / self.time_step))

        # If no seed was given in the state, set it from the seed
        if self.state[3] is None:
            object.__setattr__(self, 'state', self.state[:3] + (self.seed,))

        # Total saves is one for the start of each 'save_interval' plus one at the end of the run.
        object.__setattr__(self, 'number_of_saves', int(self.total_time / self.save_interval) + 1)

    def __str__(self) -> str:
        """Returns a human-readable, JSON-like string of all parameters."""
        # Convert the internal parameters into one dictionary
        values = asdict(self)
        # Format the dictionary and return
        return dict_to_formatted_str(values)

    @staticmethod
    def print_default_values() -> str:
        """Returns the default parameters for the Macroscale model."""
        # Create a new MacroParameters object with the default values
        default_macro_params = MacroParameters()
        # Convert to a dict, then to a formatted string, and return
        return str(default_macro_params)



