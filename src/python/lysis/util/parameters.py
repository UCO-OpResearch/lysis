"""Code for holding, storing, and reading information about a Run

This module gives a uniform way to handle the data and parameters of a given 
run. It contains classes to house these and make them accessible to the 
rest of the code. It also handles the storing and reading of parameters and 
data to/from disk.

Typical usage example:
    >>> # Create a new run
    >>> exp = Run('path/to/data')
    >>> param = {'override_parameter': 2.54, 'another_new_parameter': 32}
    >>> exp.initialize_macro_param(param)
    >>> exp.to_file()
    >>> # Load an existing run
    >>> exp = Run('path/to/data', '2022_12_27_1100')
    >>> exp.read_file()
    >>> # Access a parameter
    >>> exp.macro_params.pore_size
    >>> # Access data
    >>> exp.data.lysis_time[4][18]
"""

import inspect
import json
import logging
import os
import pkgutil
import re
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, List, Mapping, Tuple, Union

from pint import Quantity

from .constants import default_filenames, ureg, Q_
from .datastore import DataStore
from .util import dict_to_formatted_str


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2024, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.2"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# TODO: Rename this object to Run
class Run(object):
    """Houses all information about a given experimental run.

    This object contains:

    * Data location
    * Run parameters
    * Run data

    It includes methods for

    * Initializing with default parameters
    * Reading parameters from disk
    * Saving parameters to disk
    * Reading input data from disk
    * Writing result data to disk

    Args:
        data_root: The path of the folder containing datasets
        run_code: The code number of the run.
            This will be the name of the folder containing the data specific to
            this run.
            This should be a date and time in 'YYYY-MM-DD-hhmm' format
            If no code is given, one will be generated from the current date
            and time.

    Attributes:
        run_code (str): The code number of the run.
        os_path (str): The path to the folder containing this run's data
        macro_params (DataClass): A dictionary of


    Raises:
        RuntimeError: An invalid data folder was given.
    """

    def __init__(self, data_root: Union[str, bytes, os.PathLike], run_code: str = None):
        # Check if the data folder path is valid
        if not os.path.isdir(data_root):
            raise RuntimeError("Data folder not found.", data_root)
        self.os_data_root = data_root
        # If no run code was given, create a new one from the current
        # date and time.
        if run_code is None:
            self.run_code = datetime.now().strftime("%Y-%m-%d-%H%M")
        else:
            self.run_code = run_code

        # Generate the path to the run folder and the parameters file
        self.os_path = os.path.join(data_root, str(self.run_code))
        os.makedirs(self.os_path, exist_ok=True)
        self.os_param_file = os.path.join(self.os_path, "params.json")

        # TODO(bpaynter): Check if the parameters are already stored.
        #                 Don't allow parameters to be changed once stored.
        # Initialize the internal storage as empty
        # self.sequence: ExpComponent = ExpComponent.NONE
        self.micro_params = None
        self.macro_params = None
        self.data = DataStore(self.os_path, default_filenames)

    def __str__(self) -> str:
        """Gives a human-readable, formatted string of the current run's
        parameters."""
        # Convert internal storage to a dictionary
        values = self.to_dict()
        # Call the formatter and return
        return dict_to_formatted_str(values)

    def initialize_micro_param(self, params: Mapping[str, Any] = None) -> None:
        """Creates the parameters for the Microscale model.

        Parameters are set to the default values unless new values are passed
        in the params dictionary.

        This method is essentially a wrapper for the MicroParameters
        constructor.

        Args:
            params: A dictionary of parameters that differ from the default
                    values.

                For example,
                    >>> {'binding_rate': 10, 'pore_size': 3,}
        """
        if params is not None:
            self.micro_params = MicroParameters(**params)
        else:
            self.micro_params = MicroParameters()

    def initialize_macro_param(self, params: dict[str, Any] = None) -> None:
        """Creates the parameters for the Macroscale model.

        Parameters are set to the default values unless new values are passed
        in the params dictionary.

        This method is essentially a wrapper for the MacroParameters
        constructor.

        Args:
            params: A dictionary of parameters that differ from the default
                    values.

                For example,
                    >>> {'binding_rate': 10, 'pore_size': 3,}
        """
        # The macroscale model is dependent on the parameters and results of the
        # microscale model. If no microscale parameters are supplied, the macroscale
        # model cannot be initialized.
        if self.micro_params is None:
            raise RuntimeError("No Microscale parameters.")
        if params is not None:
            self.macro_params = MacroParameters(**params)
        else:
            self.macro_params = MacroParameters(micro_params=self.micro_params)

    def to_dict(self) -> dict:
        """Returns the internally stored data as a dictionary.

        Does not include system-specific information like paths.
        """
        # Initialize a dictionary of the appropriate parameters
        output = {
            "run_code": self.run_code,
            "data_filenames": None,
            "micro_params": None,
            "macro_params": None,
        }
        # Get the data filenames from the DataStore
        if self.data is not None:
            output["data_filenames"] = self.data.to_dict()
        # Convert the Microscale parameters to a dictionary
        if self.micro_params is not None:
            # Get units
            units = MacroParameters.units()
            output["micro_params"] = {}
            # Loop through the parameters
            for k, v in asdict(self.micro_params).items():
                # If the parameter is stored as a Quantity, convert it to standard units
                # and output as a string. Else, pass it as-is
                if isinstance(v, Quantity):
                    output["micro_params"][k] = str(v.to(units[k]))
                else:
                    output["micro_params"][k] = v

        # Convert the Macroscale parameters to a dictionary
        if self.macro_params is not None:
            # Get units
            units = MacroParameters.units()
            output["macro_params"] = {}
            # Loop through the parameters
            for k, v in asdict(self.macro_params).items():
                # If the parameter is stored as a Quantity, convert it to standard units
                # and output as a string. Else, pass it as-is
                if isinstance(v, Quantity):
                    output["macro_params"][k] = str(v.to(units[k]))
                else:
                    output["macro_params"][k] = v

        return output

    def to_file(self) -> None:
        """Stores the run parameters to disk.

        Creates or overwrites the params.json file in the run's data
        folder. This file will contain the current run parameters
        (including any micro- and macroscale parameters) in JSON format.
        """
        with open(self.os_param_file, "w") as file:
            # Convert the internal parameters to a dictionary and then use the
            # JSON module to save to disk.
            json.dump(self.to_dict(), file)

    def read_file(self) -> None:
        """Load the run parameters from disk.

        Raises:
            RuntimeError: No parameter file is available for this run.
        """
        # Determine whether the parameter file exists for this run
        if not os.path.isfile(self.os_param_file):
            raise RuntimeError("Run parameter file not found.")
        # Open the file
        with open(self.os_param_file, "r") as file:
            # Use the JSON library to read in the parameters as a dictionary
            params = json.load(file)
            # Initialize a datastore
            data_filenames = params.pop("data_filenames", None)
            if data_filenames is not None:
                self.data = DataStore(self.os_path, data_filenames)

            # Remove the Microscale parameters from the dictionary (if it
            # exists) and create a new MicroParameters object using its values
            micro_params = params.pop("micro_params", None)
            if micro_params is not None:
                # We are checking here to make sure that saved, dependent
                # parameters don't get passed to the MicroParameters
                # constructor
                out_micro_params = {}
                units = MicroParameters.units()
                # Find the parameters needed to initialize a new
                # MicroParameters object
                sig = inspect.signature(MicroParameters)
                # Get the keys we read from the JSON
                for k, v in micro_params.items():
                    # If that key is not needed, then toss it
                    if k not in sig.parameters:
                        continue
                    # If the parameter has units, parse it with Pint
                    if k in units:
                        out_micro_params[k] = Q_(v)
                    else:
                        out_micro_params[k] = v
                # Now unpack whatever is left in the dict and pass it to the
                # constructor
                self.micro_params = MicroParameters(**out_micro_params)
            else:
                # If there were no microscale parameters in the file, then
                # we raise a warning.
                warnings.warn(
                    "Run parameter file does not contain Microscale parameters. "
                    "Using defaults.",
                    RuntimeWarning,
                )
                self.micro_params = MicroParameters()

            # Remove the Macroscale parameters from the dictionary (if it
            # exists) and create a new MacroParameters object using its values
            macro_params = params.pop("macro_params", None)
            if macro_params is not None:
                # We are checking here to make sure that saved, dependent
                # parameters don't get passed to the MacroParameters
                # constructor

                out_macro_params = {}
                units = MacroParameters.units()
                # Find the parameters needed to initialize a new
                # MicroParameters object
                sig = inspect.signature(MacroParameters)
                # Get the keys we read from the JSON
                for k, v in macro_params.items():
                    # If that key is not needed, then toss it
                    if k not in sig.parameters:
                        continue
                    # If the parameter has units, parse it with Pint
                    if k in units:
                        out_macro_params[k] = Q_(v)
                    else:
                        out_macro_params[k] = v
                # Now unpack whatever is left in the dict and pass it to the
                # constructor
                self.macro_params = MacroParameters(**out_macro_params)
            else:
                # If there were no parameters in the file, then we leave the
                # object null.
                self.macro_params = None


################################
###  NOTE:
###  In the following dataclass definitions, do not use double-quotes (")
###  Except for the docstrings below each definition.
###  Otherwise it will mess up the units() and fortran() regexes
###
### TODO: Fix the regexes so that this is no longer a problem
################################


@dataclass(frozen=True)
class MicroParameters:
    """This will contain the parameters for the Microscale model.

    Parameters can be accessed as attributes.
    Independent parameters should only be set at initialization.
    Dependent parameters should never be set manually, but are automatically
    calculated by internal code.

    Should only be used inside a Run object.


    Example:
        >>> # Initialize using the default values
        >>> micro_params_default = MicroParameters()
        >>> # Get parameter value
        >>> micro_params_default.fibrinogen_length
        45 nanometers
        >>> # Initialize overriding some default values
        >>> p = {'fibrinogen_radius': "10 nm", 'fiber_radius': "3 um"}
        >>> micro_params_override = MicroParameters(**p)
    """

    # TODO(bpaynter): Potential future features:
    #                   * The inclusion of standard sets of Microscale
    #                       parameters (i.e., Q2, CaseA-D, etc.)

    #####################################
    # Physical Parameters
    #####################################

    fibrinogen_length: Quantity = Q_("45 nanometers")
    """The length of a fibronogen molecule.
    
    :Units: microns
    :Fortran: None"""

    # NOTE: Currently set to 1.2 nanometers to match the legacy 2.4nm protofibril radius.
    # This should be changed to 2.5 nanometers once verification is complete to match
    # Yeromonahos, 2010 doi: 10.1016/j.bpj.2010.04.059
    fibrinogen_radius: Quantity = Q_("1.2 nanometers")
    """The radius of a fibronogen molecule.
    
    :Units: microns
    :Fortran: None"""

    fiber_radius: Quantity = Q_("72.7/2 nanometers")
    """The radius of each fiber in the model.
    
    :Units: microns
    :Fortran: radius"""

    # This was changed from 2.4 nm on 2024-01-17 to match physiological values
    # Yeromonahos, 2010 doi: 10.1016/j.bpj.2010.04.059
    protofibril_radius: Quantity = field(init=False)
    """The radius of a protofibril.
    
    :Units: microns
    :Fortran: None"""

    diss_const_tPA_wPLG: Quantity = Q_("0.02 micromolar")
    """The dissociation constant of tPA, :math:`k^D_\\text{tPA}`, to fibrin 
    in the presence of PLG.
    
    :Units: micromolar
    :Fortran: KdtPAyesplg"""

    diss_const_tPA_woPLG: Quantity = Q_("0.36 micromolar")
    """The dissociation constant of tPA, :math:`k^D_\\text{tPA}`, to fibrin
    in the absence of PLG.

    :Units: micromolar
    :Fortran: KdtPAnoplg"""

    diss_const_PLG_intact: Quantity = Q_("38 micromolar")
    """The dissociation constant of PLG, :math:`k^D_\\text{PLG}`, to intact fibrin.

    :Units: micromolar
    :Fortran: KdPLGintact"""

    diss_const_PLG_nicked: Quantity = Q_("2.2 micromolar")
    """The dissociation constant of PLG, :math:`k^D_\\text{PLG}`, to nicked fibrin.

    :Units: micromolar
    :Fortran: KdPLGnicked"""

    bind_rate_tPA: Quantity = Q_("0.1 (micromolar*sec)^-1")
    """The binding rate of tPA, :math:`k^\\text{on}_\\text{tPA}`, to fibrin.

    :Units: (micromolar*sec)^-1
    :Fortran: ktPAon"""

    bind_rate_PLG: Quantity = Q_("0.1 (micromolar*sec)^-1")
    """The binding rate of PLG, :math:`k^\\text{on}_\\text{PLG}`, to fibrin.

    :Units: (micromolar*sec)^-1
    :Fortran: kPLGon"""

    conc_free_PLG: Quantity = Q_("2 micromolar")
    """The concentration of free plasminogen.
    
    :Units: micromolar
    :Fortran: freeplg"""

    deg_rate_fibrin: Quantity = Q_("5 sec^-1")
    """The plasmin-mediated rate of fibrin degradation.
    
    :Units: sec^-1
    :Fortran: kdeg"""

    unbind_rate_PLG_intact: Quantity = field(init=False)
    """The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, 
    from intact fibrin.

    :Units: sec^-1
    :Fortran: kplgoff"""

    unbind_rate_PLG_nicked: Quantity = field(init=False)
    """The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, 
    from nicked fibrin.

    :Units: sec^-1
    :Fortran: kplgoffnick"""

    unbind_rate_PLi: Quantity = Q_("57.6 sec^-1")
    """The unbinding rate of PLi, :math:`k^\\text{off}_\\text{PLi}`, 
    from fibrin.

    :Units: sec^-1
    :Fortran: kplioff"""

    unbind_rate_tPA_wPLG: Quantity = field(init=False)
    """The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, 
    from fibrin in the presence of PLG.

    :Units: sec^-1
    :Fortran: kaoff12"""

    unbind_rate_tPA_woPLG: Quantity = field(init=False)
    """The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, 
    from fibrin in the absence of PLG.

    :Units: sec^-1
    :Fortran: kaoff10"""

    activation_rate_PLG: Quantity = Q_("0.1 sec^-1")
    """The catalytic rate constant, :math:`k_\\text{cat}^\\text{ap}`, 
    for activation of PLG into PLI.
    
    :Units: sec^-1
    :Fortran: kapcat"""

    exposure_rate_binding_site: Quantity = Q_("5 sec^-1")
    """The catalytic rate constant, :math:`k_\\text{cat}^\\text{n}`, 
    for the PLi-mediated rate of exposure of new binding sites.

    :Units: sec^-1
    :Fortran: kncat"""

    protein_per_fiber: Quantity = field(init=False)
    """The fraction of protein in each fiber (by volume?)
    
    :Units: %
    :Fortran: None"""

    fibrin_conc_per_fiber: Quantity = field(init=False)
    """The concentration of fibrin in each fiber

    :Units: micromolar
    :Fortran: None"""

    binding_sites: Quantity = field(init=False)  # int = 427
    """Concentration of binding sites.
     
    :Units: micromolar
    :Fortran: bs"""

    #####################################
    # Model Parameters
    #####################################

    nodes_in_row: int = 7
    """The number of protofibrils in one row of the lattice inside one
    fiber.
    
    :Units: None
    :Fortran: nodes"""

    snap_proportion: float = 2.0 / 3.0
    """The proportion of doublets that need to be degraded before the
    fiber snaps.
    
    :Units: None
    :Fortran: snap_proportion"""

    #####################################
    # Mechanism Parameters
    #####################################

    simulations: int = 50_000
    """The number of independent trials run in the microscale model.
    
    :Units: None
    :Fortran: simulations"""

    seed: int = 0
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    #####################################
    # Data Parameters
    #####################################

    output_data: List[str] = field(init=False)
    """The data output by the Microscale model."""

    #####################################
    # Code Parameters
    #####################################

    micro_version: str = "micro_rates"
    """A string identifying which version of the Microscale model is being run."""

    log_lvl: int = logging.WARNING
    """How much debugging information to write out to the console

    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        """This method calculates the dependent parameters once the
        MicroParameters object is created. It is automatically called by the
        DataClass.__init__()"""

        # These names must be elements of the Run's DataStore
        object.__setattr__(
            self,
            "output_data",
            [
                "lysis_complete_time",  # Fortran: lysis_time
                "tPA_leaving_time",  # Fortran: tPA_time
                "PLi_generated",  # Fortran: Plasmin
                "lysis_completed",  # Fortran: lysiscomplete
                "tPA_kinetic_unbound",  # Fortran: tPAunbind
                "tPA_forced_unbound",  # Fortran: tPAPLiunbd
                "tPA_still_bound",  # Fortran: ltPA
                "first_PLi",  # Fortran: firstPLi
            ],
        )

        # One protofibril is two fibrinogens
        object.__setattr__(
            self,
            "protofibril_radius",
            2 * self.fibrinogen_radius,
        )

        # The dissociation constant is the unbinding rate over the binding rate
        object.__setattr__(
            self,
            "unbind_rate_PLG_intact",
            self.bind_rate_PLG * self.diss_const_PLG_intact,
        )

        # The dissociation constant is the unbinding rate over the binding rate
        object.__setattr__(
            self,
            "unbind_rate_PLG_nicked",
            self.bind_rate_PLG * self.diss_const_PLG_nicked,
        )

        # The dissociation constant is the unbinding rate over the binding rate
        object.__setattr__(
            self,
            "unbind_rate_tPA_wPLG",
            self.bind_rate_tPA * self.diss_const_tPA_wPLG,
        )

        # The dissociation constant is the unbinding rate over the binding rate
        object.__setattr__(
            self,
            "unbind_rate_tPA_woPLG",
            self.bind_rate_tPA * self.diss_const_tPA_woPLG,
        )

        # The fraction of fiber which is protein
        # Equation on page S2 from Bannish, et. al. 2017
        # https://doi.org/10.1038/s41598-017-06383-w
        object.__setattr__(
            self,
            "protein_per_fiber",
            (
                self.nodes_in_row**2
                / (self.fibrinogen_length / 2 * ureg.pi * self.fiber_radius**2)
                * (self.fibrinogen_length / 2)
                * (ureg.pi * self.protofibril_radius**2)
            ).to("%"),
        )

        # The fibrin concentration of each fiber
        # Equation on page S2 from Bannish, et. al. 2017
        # https://doi.org/10.1038/s41598-017-06383-w
        object.__setattr__(
            self,
            "fibrin_conc_per_fiber",
            (
                self.nodes_in_row**2
                / (self.fibrinogen_length / 2 * ureg.pi * self.fiber_radius**2)
                / ureg.avogadro_constant
            ).to("micromolar"),
        )

        # Calculate the concentration of binding sites per fiber
        # Calculation on page S3 from Bannish, et. al. 2017
        # https://doi.org/10.1038/s41598-017-06383-w
        object.__setattr__(
            self,
            "binding_sites",
            4
            * (self.nodes_in_row - 1)
            / self.nodes_in_row**2
            * self.fibrin_conc_per_fiber,
        )

    @staticmethod
    def units():
        """Returns a dictionary whose keys are the names of all parameters (both micro- and macroscale)
        that have units. The values in the dictionary are those units.
        These values are parsed from the docstrings in this file.
        """
        # Get the text of this source code file
        text = pkgutil.get_data(__name__, "parameters.py")
        # Initialize the regex for finding :Units: tags in the docstrings
        pattern = re.compile(
            r"[\r\n]"  # The start of a line
            + r"^\s{4}"  # four spaces (one tab)
            + r"([a-zA-Z0-9_]+)"  # First capture group, the name of the parameter
            + r":.*[\r\n]*"  # the rest of that line
            + r".*\"\"\""  # a triple-quote
            + r"[^\"]*"  # text that is not a douple-quote
            + r":Units:\s"  # the units tag
            + r"([^\n\r]+)"  # Second capture group, the units
            + r"[\r\n]",  # A line break
            re.M,
        )
        units = {}
        # Find all matches
        matches = re.findall(pattern, text.decode("utf-8"))

        # Check that the units are not none.
        for match in matches:
            if match[1] != "None":
                units[match[0]] = match[1]

        return units

    @staticmethod
    def fortran_names():
        """Returns a dictionary whose keys are the names of all parameters (both micro- and macroscale)
        that have equivalents in Fortran. The values in the dictionary are the names of the equivalent
        Fortran variable names.
        These values are parsed from the docstrings in this file.
        """
        # Get the text of this source code file
        text = pkgutil.get_data(__name__, "parameters.py")
        # Initialize the regex for finding :Units: tags in the docstrings
        pattern = re.compile(
            r"[\r\n]+"  # The start of a line
            + r"^\s{4}"  # four spaces (one tab)
            + r"([a-zA-Z0-9_]+)"  # First capture group, the name of the parameter
            + r":.*[\r\n]*"  # the rest of that line
            + r".*\"\"\""  # a triple-quote
            + r"[^\"]*"  # text that is not a douple-quote
            + r":Fortran:\s"  # the Fortran tag
            + r"([\w_]+(-1)?)"  # Second capture group, the fortran name, possibly with a -1
            + r"(\s=[^\"]*)?"  # Third capture group (optional), a formula for the parameter in Fortran
            + r"\"\"\"",  # A line break
            re.M,
        )
        names = {}
        # Find all matches
        matches = re.findall(pattern, text.decode("utf-8"))
        # Check that the units are not none.
        for match in matches:
            if match[1] != "None":
                names[match[0]] = match[1]

        return names

    def __str__(self) -> str:
        """Returns a human-readable, JSON-like string of all parameters."""
        # Convert the internal parameters into one dictionary
        values = asdict(self)
        # Format the dictionary and return
        return dict_to_formatted_str(values)

    @staticmethod
    def print_default_values() -> str:
        """Returns the default parameters for the Microscale model."""
        # Create a new MicroParameters object with the default values
        default_micro_params = MicroParameters()
        # Convert to a dict, then to a formatted string, and return
        return str(default_micro_params)


@dataclass(frozen=True)
class MacroParameters:
    """Contains parameters for the Macroscale model.

    Parameters can be accessed as attributes.
    Independent parameters should only be set at initialization.
    Dependent parameters should never be set manually, but are automatically
    calculated by internal code.

    Should only be used inside a Run object.


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

    pore_size: Quantity = Q_("1.0135 um")
    """Pore size (distance between fibers/nodes)
    
    :Units: centimeters
    :Fortran: delx"""

    diffusion_coeff: Quantity = Q_("5.0e-7 cm^2/s")
    """Diffusion coefficient
    
    :Units: cm^2/s
    :Fortran: Diff"""

    # TODO(bpaynter): This value should derive from MicroParameters
    forced_unbind: float = 0.0852
    """Fraction of times tPA was forced to unbind in microscale model.
    
    :Units: None
    :Fortran: frac_forced"""

    # TODO(bpaynter): This value should derive from MicroParameters
    # TODO(bpaynter): Rename to average_bound_time
    average_bound_time: Quantity = Q_("27.8 sec")
    """This is the average time a tPA molecule stays bound to fibrin. 
    For now I'm using 27.8 to be 1/0.036, the value in the absence of PLG.
    
    :Units: seconds
    :Fortran: avgwait = 1/koff"""

    #####################################
    # Model Parameters
    #####################################

    cols: int = 93
    """The number of lattice nodes in each (horizontal) row
    
    :Units: None
    :Fortran: N"""

    # TODO(bpaynter): 'rows' and 'fiber_rows' should be switched so that
    #                 'fiber_rows' is the independent variable.
    rows: int = 121
    """The number of lattice nodes in each (vertical) column
    
    :Units: None
    :Fortran: F"""

    fiber_rows: int = field(init=False)
    """The number of rows containing fibrin
    
    :Units: None
    :Fortran: Fhat"""

    empty_rows: int = 29 - 1
    """The number of fibrin-free rows at the top of the grid.
    
    Equivalent to 'first_fiber_row', which is the 1st node in vertical 
    direction containing fibers.
    So if first_fiber_row = 10, then rows 0-9 have no fibers, there's one more 
    row of fiber-free planar vertical edges, and then the row with index 
    'first_fiber_row' (e.g. 11th) is a full row of fibers.
    
    
    :Units: None
    :Fortran: Ffree-1"""

    empty_edges: int = field(init=False)
    """The number of edges without fibrin.
    Also the 1-D index of the last edge without fibrin when 1-indexing
    
    This is probably unnecessary when using a 2-D data structure, but is kept 
    for historical reasons.
    
    
    :Units: None
    :Fortran: enoFB"""

    full_row: int = field(init=False)
    """Edges in a full row of nodes
    
    :Units: None
    :Fortran: None"""

    xz_row: int = field(init=False)
    """Number of all x- and z-edges in a row
    
    :Units: None
    :Fortran: None"""

    total_edges: int = field(init=False)
    """The total number of edges in the model
    
    :Units: None
    :Fortran: num"""

    total_fibers: int = field(init=False)
    """The total number of fibers in the model

    :Units: None
    :Fortran: None"""

    total_molecules: int = 43074
    """The total number of tPA molecules:
    
        * 43074 is Colin's [tPA]=0.6 nM
        * 86148 is Colin's [tPA]=1.2 nM
        
    :Units: None
    :Fortran: M"""

    moving_probability: float = 0.2
    """The probability of moving.
    
    Make sure it is small enough that we've converged.
    
    :Units: None
    :Fortran: q"""

    #####################################
    # Mechanism Parameters
    #####################################

    simulations: int = 10
    """The number of independent simulations to be run
    
    :Units: None
    :Fortran: simulations"""

    total_time: Quantity = Q_("20 min")
    """Total running time for model.
     
    :Units: seconds
    :Fortran: tf"""

    time_step: float = field(init=False)
    """The length of one timestep.
    
    :Units: seconds
    :Fortran: tstep"""

    total_time_steps: int = field(init=False)
    """The total number of timesteps.
    
    :Units: None
    :Fortran: num_t"""

    seed: int = 0  # -2137354075
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    state: Tuple[int, int, int, int] = field(init=False)
    """State for the random number generator.
    
    :Units: None
    :Fortran: state"""

    #####################################
    # Data Parameters
    #####################################

    input_data: List[str] = field(init=False)
    """The data (from the Microscale model) required to run the Macroscale 
    model."""

    output_data: List[str] = field(init=False)
    """The data output by the Macroscale model."""

    save_interval: Quantity = Q_("10 sec")
    """How often to record data from the model.
    
    :Units: sec
    :Fortran: save_interval"""

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
    This will impact performance negatively.
    This currently does nothing.

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
        # These names must be elements of the Run's DataStore
        object.__setattr__(
            self,
            "input_data",
            [
                "unbinding_time",  # Fortran: tsec1
                # 'leaving_time',           # Fortran: CDFtPA
                "lysis_time_dist",  # Fortran: lysismat
                "total_lyses",  # Fortran: lenlysismat
            ],
        )
        # These names must be elements of the Run's DataStore
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

        # The number of fiber rows is the total rows minus the empty ones
        object.__setattr__(self, "fiber_rows", self.rows - self.empty_rows)

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
        # number of edges in the fibrin-free region.
        # The total rows in the fibrin-free region is equal to the (0-based)
        # index of the first fiber row.
        # The total edges in this region is one full row of edges for each row
        object.__setattr__(self, "empty_edges", self.full_row * self.empty_rows)

        # Equation (2.4) page 25 from Bannish, et. al. 2014
        # https://doi.org/10.1093/imammb/dqs029
        object.__setattr__(
            self,
            "time_step",
            (
                self.moving_probability
                * self.pore_size**2
                / (12 * self.diffusion_coeff)
            ).to_reduced_units(),
        )

        # Total timesteps is total time divided by length of one timestep
        object.__setattr__(
            self, "total_time_steps", int(self.total_time / self.time_step)
        )

        # Set the state
        object.__setattr__(self, "state", (129281, 362436069, 123456789, self.seed))

        # Total saves is one for the start of each 'save_interval' plus one at
        # the end of the run.
        object.__setattr__(
            self, "number_of_saves", int(self.total_time / self.save_interval) + 1
        )

    @staticmethod
    def units():
        """Returns a dictionary whose keys are the names of all parameters (both micro- and macroscale)
        that have units. The values in the dictionary are those units.
        These values are parsed from the docstrings in this file.
        """
        # Get the text of this source code file
        text = pkgutil.get_data(__name__, "parameters.py")
        # Initialize the regex for finding :Units: tags in the docstrings
        pattern = re.compile(
            r"[\r\n]"  # The start of a line
            + r"^\s{4}"  # four spaces (one tab)
            + r"([a-zA-Z0-9_]+)"  # First capture group, the name of the parameter
            + r":.*[\r\n]*"  # the rest of that line
            + r".*\"\"\""  # a triple-quote
            + r"[^\"]*"  # text that is not a douple-quote
            + r":Units:\s"  # the units tag
            + r"([^\n\r]+)"  # Second capture group, the units
            + r"[\r\n]",  # A line break
            re.M,
        )
        units = {}
        # Find all matches
        matches = re.findall(pattern, text.decode("utf-8"))

        # Check that the units are not none.
        for match in matches:
            if match[1] != "None":
                units[match[0]] = match[1]

        return units

    @staticmethod
    def fortran_names():
        """Returns a dictionary whose keys are the names of all parameters (both micro- and macroscale)
        that have equivalents in Fortran. The values in the dictionary are the names of the equivalent
        Fortran variable names.
        These values are parsed from the docstrings in this file.
        """
        # Get the text of this source code file
        text = pkgutil.get_data(__name__, "parameters.py")
        # Initialize the regex for finding :Units: tags in the docstrings
        pattern = re.compile(
            r"[\r\n]+"  # The start of a line
            + r"^\s{4}"  # four spaces (one tab)
            + r"([a-zA-Z0-9_]+)"  # First capture group, the name of the parameter
            + r":.*[\r\n]*"  # the rest of that line
            + r".*\"\"\""  # a triple-quote
            + r"[^\"]*"  # text that is not a douple-quote
            + r":Fortran:\s"  # the Fortran tag
            + r"([\w_]+(-1)?)"  # Second capture group, the fortran name, possibly with a -1
            + r"(\s=[^\"]*)?"  # Third capture group (optional), a formula for the parameter in Fortran
            + r"\"\"\"",  # A line break
            re.M,
        )
        names = {}
        # Find all matches
        matches = re.findall(pattern, text.decode("utf-8"))
        # Check that the units are not none.
        for match in matches:
            if match[1] != "None":
                names[match[0]] = match[1]

        return names

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
