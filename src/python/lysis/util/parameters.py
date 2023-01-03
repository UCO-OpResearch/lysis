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
    >>> exp = Experiment('path/to/data', '2022-12-27-1100')
    >>> exp.read_file()
    >>> # Access a parameter
    >>> exp.macro_params['pore_size']

    >>> # Read data from disk
    >>> exp.read_macro_input_data()
    >>> # Access data
    >>> exp.data.lysis_time[4][18]          # type: ignore
"""

import errno
import json
import logging
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, AnyStr, Mapping, Tuple

import numpy as np

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
        micro_params (dict): A dictionary of


    Raises:
        RuntimeError: An invalid data folder was given.
    """
    def __init__(self, data_root: AnyStr, experiment_code: str = None):
        # Check if the data folder path is valid
        if not os.path.isdir(data_root):
            raise RuntimeError('Data folder not found.')
        self.os_data_root = data_root
        # If no experiment code was given, create a new one from the current date and time.
        if experiment_code is None:
            now = datetime.now()
            self.experiment_code = str.join('', [str(now.year),
                                                 '_',
                                                 str(now.month),
                                                 '_',
                                                 str(now.day),
                                                 '_',
                                                 str(now.hour),
                                                 str(now.day)
                                                 ])
        else:
            self.experiment_code = experiment_code

        # Generate the path to the experiment folder and the parameters file
        self.os_path = os.path.join(data_root, str(self.experiment_code))
        self.os_param_file = os.path.join(self.os_path, 'params.json')

        # Initialize the internal storage as empty
        self.micro_params: Mapping[str, Any] | None = None
        self.macro_params = None
        self.data = Data()

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

    def to_dict(self) -> dict:
        """Returns the internally stored data as a dictionary.

        Does not include system-specific information like paths.
        """
        # Initialize a dictionary of the appropriate parameters
        output = {'experiment_code': self.experiment_code,
                  'micro_params': None,
                  'macro_params': None
                  }
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

    # TODO(bpaynter): This is broken since the conversion to DataClass and needs to be fixed.
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
            # Remove the Macroscale parameters from the dictionary (if it exists)
            # and create a new MacroParameters object using its values
            macro_params = params.pop('macro_params', None)
            if macro_params is not None:
                self.macro_params = MacroParameters(**macro_params)
            else:
                self.macro_params = None

    # TODO(bpaynter): Redo with a DataClass
    def read_macro_input_data(self) -> None:
        """Reads matrices used by the macroscale model from disk.

        Performs some cleanup on data files if they are in text format:

        * Will convert total_lyses to integer
        * Will transpose the lysis_time matrix
        * Will trim the lysis_time matrix to remove excess null (6000) values

        This cleanup should be unnecessary once the whole model (incl microscale) is converted to Python and NumPy.

        Raises:
            RuntimeError: Macroscale parameters are empty.
            FileNotFoundError: A matrix data file is not available.
            NotImplementedError: Support for reading NumPy files needs to be added.
            AttributeError: A non-supported file type is requested.
        """
        # Check if the macroscale parameters are available
        if self.macro_params is None:
            raise RuntimeError('Macro Parameters not initialized.')
        # Get the requested matrices from the parameters
        for array, filename in self.macro_params.data_files.items():
            # Check if the requested filename exists
            if not os.path.isfile(os.path.join(self.os_path, filename)):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.join(self.os_path, filename))
            # If the file is stored in text format (with the extension '.dat')
            # This is the format used by the original Fortran/Matlab code
            if os.path.splitext(filename)[1] == '.dat':
                # The total_lyses matrix needs to be converted to integer
                if array == 'total_lyses':
                    self.data[array] = np.loadtxt(os.path.join(self.os_path, filename), dtype=int, converters=float)
                # The lysis_time matrix needs to be transposed
                elif array == 'lysis_time':
                    self.data[array] = np.loadtxt(os.path.join(self.os_path, filename)).T
                # Else just read in the array normally
                else:
                    self.data[array] = np.loadtxt(os.path.join(self.os_path, filename))
            # If the file is stored in NumPy format (with the extension .npy)
            # This will be the format used by the Python model once complete
            # TODO(bpaynter): Implement, consistent with the rest of the model
            elif os.path.splitext(filename)[1] == '.npy':
                raise NotImplementedError('Support for NumPy files is yet to be implemented.')
            # Other file types/extensions are not supported (or planned to be) at this time
            else:
                raise AttributeError('Non-supported file type.')
        # The lysis_time matrix is currently stored with a full 500 columns.
        # The total_lyses matrix contains the number of finite values in each row of the lysis_time matrix
        # At the end of finite data, the rest of the row is filled with the value 6000 to indicate null or infinite.
        # Thus, we trim the lysis_time matrix to the length of the longest row of finite data (plus one).
        self.data.lysis_time = np.delete(self.data.lysis_time, np.s_[self.data.total_lyses.max():], 1)  # type: ignore

    def save_data(self, data: Mapping[AnyStr, np.ndarray], as_text: bool = False) -> None:
        """Save a selection of NumPy matrices to disk.

        A selection of NumPy arrays is passed in a dictionary (keyed with the name of the matrix/file).
        Depending on the value of the optional as_text flag, these are stored in separate files
        as either text (with the extension .dat) or in NumPy format (with the extension .npy).

        Args:
            data: A mapping of matrix/file names to the NumPy arrays themselves
            as_text: True if the matrices are to be saved in text format, False if they are to be saved in NumPy format.
        """
        # TODO(bpaynter): It needs to be checked whether the text files output by this method can be read by
        #   the existing Fortran and Matlab code.
        for name, array in data.items():
            if as_text:
                # Save as a text file with the extension .dat
                np.savetxt(os.path.join(self.os_path, name + '.dat'),
                           array,   # type: ignore  # Not really sure. Issue with the Typing in NumPy?
                           newline=os.linesep)
            else:
                # Save as a NumPy file with the extension .npy
                np.save(os.path.join(self.os_path, name))


# TODO(bpaynter): Needs to be implemented. This implementation should include:
#                   * Standard Microscale parameters
#                   * The inclusion of standard sets of Microscale parameters (i.e., Q2, CaseA-D, etc.)
#                   * The modification of the Microscale Fortran code to read/write JSON parameters
#                   * The modification of the MacroParameters class to derive the appropriate parameters
#                       from the MicroParameters class
class MicroParameters(dict):
    """This will contain the parameters for the Microscale model.
    This will need to be implemented and the Fortran microscale code modified to work with it.

    Once implemented, many Macroscale parameters will need to be redefined so that they derive from the appropriate
    Microscale parameters.
    """
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

    rows: int = 93
    """The number of lattice nodes in each (horizontal) row
    
    :Units: nodes
    :Fortran: N"""

    cols: int = 121
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

    first_fiber_row: int = 29 - 1
    """1st node in vertical direction containing fibers.
    So if first_fiber_row = 10, then rows 0-9 have no fibers, there's one more row of fiber-free planar vertical edges,
    and then the row with index 'first_fiber_row' (e.g. 11th) is a full row of fibers
    
    :Units: nodes
    :Fortran: Ffree-1"""

    last_ghost_edge: int = field(init=False)
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

    total_trials: int = 10
    """The number of independent trials to be run
    
    :Units: trials
    :Fortran: stats"""

    total_time: int = 10 * 60
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

    data_files: dict[str, str] = None
    """Data file names"""

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

    duplicate_fortran: bool = True
    """Whether the Python code should follow the Fortran code step-by-step.
    Theoretically, with this set to "True", both sets of code will produce the exact same output.
    This will impact performance negatively.

    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        if self.data_files is None:
            object.__setattr__(self, 'data_files', {
                                                        'unbinding_time': "tsectPA.dat",        # Fortran: tsec1
                                                        # 'leaving_time': "tPAleave.dat",       # Fortran: CDFtPA
                                                        'lysis_time':     "lysismat.dat",       # Fortran: lysismat
                                                        'total_lyses':    "lenlysisvect.dat",   # Fortran: lenlysismat
                                                    })
        # A full row of the fiber grid contains a 'right', 'up', and 'out' edge for each node,
        # except the last node which contains no 'right' edge.
        object.__setattr__(self, 'full_row', 3 * self.cols - 1)

        # A full row of 'right' and 'out' edges is two per node, except the last node which has no 'right' edge.
        object.__setattr__(self, 'xz_row', 2 * self.cols - 1)

        # The total number of edges in the grid is a full_row for each, except the last row which has no 'up' edges.
        object.__setattr__(self, 'total_edges', self.full_row * (self.rows - 1) + self.xz_row)

        # The 1-D index of the last edge in the fibrin-free region is the total number of edges in the
        # fibrin-free region -1.
        # The total rows in the fibrin-free region is equal to the (0-based) index of the first fiber row
        # The total edges in this region is one full row of edges for each row
        object.__setattr__(self, 'last_ghost_edge', self.full_row * self.first_fiber_row - 1)

        # Equation (2.4) page 25 from Bannish, et. al. 2014
        # https://doi.org/10.1093/imammb/dqs029
        object.__setattr__(self, 'time_step', (self.moving_probability
                                               * self.pore_size ** 2
                                               / (12 * self.diffusion_coeff)))

        # Total timesteps is total time divided by length of one timestep
        object.__setattr__(self, 'total_time_steps', self.total_time / self.time_step)

        # If no seed was given in the state, set it from the seed
        if self.state[3] is None:
            object.__setattr__(self, 'state', self.state[:3] + (self.seed,))

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


class Data(dict):
    """Encapsulates a dictionary so that its elements can be accessed as properties.

    e.g.,
        >>> my_data = Data()
        >>> my_data['test'] = 'hat'
        >>> my_data.test                        # type: ignore # Python allows this, but is a little confused
        'hat'
    """
    def __init__(self):
        """Creates a blank Data object."""
        # Call the super-constructor
        super(Data, self).__init__()
        # This is where the wizardry happens.
        # We define the classes parameters as the dictionary itself
        self.__dict__ = self


def dict_to_formatted_str(d: Mapping[AnyStr, Any]) -> str:
    """Converts a dictionary into a formatted, JSON-like string.

    Align keys and values, including the alignment of sub-dicts.

    e.g.,
        >>> measurements = {'Hat': 'Large', 'Pants': {'Waist': 40, 'Inseam': 38}}
        >>> dict_to_formatted_str(measurements)
        Hat   : Large
        Pants : Waist  : 40
                Inseam : 38

    Args:
        d (dict): A dictionary with string-like keys.
    """
    # Initialize the output string
    output = ''
    # Get the system line separator
    nl = os.linesep
    # Determine the longest key
    key_len = 0
    for k in d.keys():
        key_len = max(len(k), key_len)
    # Add a space to the key length for clearance
    key_len += 1
    # Determine the tab space for sub-dictionaries
    tab = ' ' * (key_len + 2)
    # Iterate over the dictionary
    for k, v in d.items():
        # If the value is a dictionary itself...
        if isinstance(v, dict):
            # Then we write the key
            output += f'{k:<{key_len}}: '
            # Then call this function recursively on the value, and indent it appropriately.
            output += tab.join(dict_to_formatted_str(v).splitlines(True))
        # Otherwise we just print the key and its value.
        else:
            output += f'{k:<{key_len}}: {v}' + nl
    return output
