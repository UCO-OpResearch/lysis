"""Parameter definitions for microscale and macroscale simulations.

This module defines the parameter classes used throughout the lysis simulation
package. It provides:

- Base Parameters class with serialization and metadata extraction
- MicroParameters: Configuration for microscale (fiber-level) simulation
- MacroParameters: Configuration for macroscale (clot-level) simulation

The parameter classes use frozen dataclasses to ensure immutability after
initialization. Some parameters are independent (set by user), while others
are dependent (calculated automatically in __post_init__).

Parameter Metadata System
--------------------------

The module includes a regex-based system for extracting parameter metadata from
docstrings. Each parameter can include special tags:

- ``:Units: <unit_string>`` - Physical units for the parameter (parsed by units())
- ``:Fortran: <fortran_name>`` - Equivalent Fortran variable name (parsed by fortran_names())

These tags are parsed at runtime using regular expressions to build metadata
dictionaries. **Do not modify the format of these tags** as it will break the
parsing system.

Example parameter docstring format::

    fiber_radius: Quantity = Q_("72.7/2 nanometers")
    \"\"\"The radius of each fiber in the model.

    :Units: microns
    :Fortran: radius\"\"\"

Unit Handling
-------------

Parameters with physical units use Pint Quantity objects, which provide:
- Automatic unit conversion and validation
- Dimensional analysis
- Human-readable representation

The ureg (UnitRegistry) and Q_ (Quantity constructor) are imported from
the constants module.

Usage Example
-------------

Creating parameter sets::

    >>> # Use default values
    >>> micro = MicroParameters()
    >>> # Override specific parameters
    >>> micro_custom = MicroParameters(
    ...     fiber_radius=Q_("50 nm"),
    ...     micro_simulations=100000
    ... )
    >>> # Access parameter values
    >>> micro_custom.fiber_radius
    <Quantity(50, 'nanometer')>
    >>> # Get Fortran equivalent names
    >>> MicroParameters.fortran_names()['fiber_radius']
    'radius'
"""

import inspect
import logging
import pkgutil
import re
import warnings
from dataclasses import asdict, dataclass, field
from typing import List, Tuple, Type, TypeVar

from pint import Quantity

from .constants import ureg, Q_
from ..tools.util import dict_to_formatted_str


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.2"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


################################
###  NOTE:
###  In the following dataclass definitions, do not use double-quotes (")
###  Except for the docstrings below each definition.
###  Otherwise it will mess up the units() and fortran() regexes
###
### TODO: Fix the regexes so that this is no longer a problem
################################


@dataclass(frozen=True)
class Parameters:
    """Base class for simulation parameter sets.

    Provides common functionality for both MicroParameters and MacroParameters:
    - String representation (__str__)
    - Dictionary conversion (to_dict, to_basedict)
    - Metadata extraction from docstrings (units, fortran_names)
    - Deserialization from dictionaries (parse_from_basedict)

    This class should not be instantiated directly. Use MicroParameters or
    MacroParameters instead.

    The frozen=True argument makes all instances immutable after creation,
    ensuring parameter consistency throughout a simulation run.
    """

    def __str__(self) -> str:
        """Return a human-readable, JSON-like string representation of all parameters.

        Converts the parameter dataclass to a formatted string suitable for
        display or logging. Uses JSON-like formatting for readability.

        :return: Formatted string representation of all parameters with their values
        :rtype: str
        """
        # Convert the internal parameters into one dictionary
        values = asdict(self)
        # Format the dictionary and return
        return dict_to_formatted_str(values)

    def to_basedict(self) -> dict[str, int | float | str]:
        """Convert parameters to a dictionary with base Python types only.

        Converts all Quantity objects to their string representations (with units),
        and passes through other values as-is. Nested dictionaries are skipped.
        This is useful for serialization to JSON or other formats that don't
        support custom objects.

        :return: Dictionary mapping parameter names to base type values (int, float,
            or string representation of Quantity with units)
        :rtype: dict[str, int | float | str]
        """
        # Get units
        units = self.units()
        output = {}
        # Loop through the parameters
        for k, v in asdict(self).items():
            # If the parameter is stored as a Quantity, convert it to standard units
            # and output as a string. Else, pass it as-is
            if isinstance(v, dict):
                continue
            elif isinstance(v, Quantity):
                output[k] = str(v.to(units[k]))
            else:
                output[k] = v
        return output

    def to_dict(self) -> dict[str, int | float | str | Quantity]:
        """Convert parameters to a dictionary preserving all types.

        Returns a dictionary representation of the parameters with all original
        types preserved, including Quantity objects. This is a thin wrapper around
        dataclasses.asdict().

        :return: Dictionary mapping parameter names to their values (may include
            Quantity objects, strings, ints, floats, etc.)
        :rtype: dict[str, int | float | str | Quantity]
        """
        return asdict(self)

    @staticmethod
    def units():
        """Extract physical units for all parameters via regex parsing of docstrings.

        Parses this source file's docstrings to find all ``:Units:`` tags and
        builds a dictionary mapping parameter names to their unit strings. Only
        parameters with non-None units are included.

        The regex pattern matches parameter field docstrings in the format::

            parameter_name: type = value
            \"\"\"Description.

            :Units: unit_string
            ...\"\"\"

        :return: Dictionary mapping parameter names to unit strings (e.g.,
            {'fiber_radius': 'microns', 'total_time': 'seconds'})
        :rtype: dict[str, str]

        Warning:
            This method uses regex to parse source code. Do not modify the
            format of ``:Units:`` tags in parameter docstrings.
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
        """Extract Fortran equivalent names for parameters via regex parsing of docstrings.

        Parses this source file's docstrings to find all ``:Fortran:`` tags and
        builds a dictionary mapping Python parameter names to their Fortran
        equivalents. Only parameters with non-None Fortran names are included.

        The regex pattern matches parameter field docstrings in the format::

            parameter_name: type = value
            \"\"\"Description.

            :Fortran: fortran_name\"\"\"

        Special suffixes in Fortran names:
        - ``-1``: Indicates 0-based to 1-based index conversion needed
        - Formula after ``=``: Optional Fortran expression for the parameter

        :return: Dictionary mapping Python parameter names to Fortran variable
            names (e.g., {'fiber_radius': 'radius', 'empty_rows': 'Ffree-1'})
        :rtype: dict[str, str]

        Warning:
            This method uses regex to parse source code. Do not modify the
            format of ``:Fortran:`` tags in parameter docstrings.
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

    @classmethod
    def print_default_values(cls) -> str:
        """Return a formatted string of all default parameter values.

        Creates a new instance of the parameter class using all default values,
        then converts it to a human-readable string representation.

        :return: Formatted string showing all default parameter values
        :rtype: str
        """
        # Create a new Parameters object with the default values
        default_params = cls()
        # Convert to a dict, then to a formatted string, and return
        return str(default_params)

    T = TypeVar("T")

    @classmethod
    def parse_from_basedict(
        cls: Type[T], base_params: dict[str, int | float | str]
    ) -> T:
        """
        Creates a Parameters object from a dict of base type values.

        Any values will be checked to see if they need units and, if so,
        will be parsed into Quantity objects with Pint.

        :param base_params: A dictionary of parameter values in base types only.
        :type base_params: dict[str, int  |  float  |  str]
        :return: A dataclass of parameters
        :rtype: MicroParameters | MacroParameters
        """
        # We are checking here to make sure that saved, dependent
        # parameters don't get passed to the MicroParameters
        # constructor
        quant_params = {}
        units = cls.units()
        # Find the parameters needed to initialize a new
        # Parameters object
        sig = inspect.signature(cls)
        # Get the keys from the string dict
        for k, v in base_params.items():
            # If that key is not needed, then toss it
            if k not in sig.parameters:
                continue
            # If the parameter has units, parse it with Pint
            if k in units:
                if isinstance(v, int) or isinstance(v, float):
                    warnings.warn(
                        f"Parameter {k} has no units. Assuming {units[k]}.",
                        RuntimeWarning,
                    )
                    quant_params[k] = Q_(v, units[k])
                else:
                    quant_params[k] = Q_(v)
            else:
                quant_params[k] = v
        return cls(**quant_params)


@dataclass(frozen=True)
class MicroParameters(Parameters):
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

    Calculated as 2 × fibrinogen_radius (one protofibril = two fibrinogens).
    This is a dependent parameter computed in __post_init__().

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
    :Fortran: kplgon"""

    conc_free_PLG: Quantity = Q_("2 micromolar")
    """The concentration of free plasminogen.
    
    :Units: micromolar
    :Fortran: freeplg"""

    deg_rate_fibrin: Quantity = Q_("5 sec^-1")
    """The plasmin-mediated rate of fibrin degradation.
    
    :Units: sec^-1
    :Fortran: kdeg"""

    unbind_rate_PLG_intact: Quantity = field(init=False)
    """The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, from intact fibrin.

    Calculated from dissociation constant and binding rate. This is a dependent
    parameter computed in __post_init__().

    :Units: sec^-1
    :Fortran: kplgoff"""

    unbind_rate_PLG_nicked: Quantity = field(init=False)
    """The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, from nicked fibrin.

    Calculated from dissociation constant and binding rate. This is a dependent
    parameter computed in __post_init__().

    :Units: sec^-1
    :Fortran: kplgoffnick"""

    unbind_rate_PLi: Quantity = Q_("57.6 sec^-1")
    """The unbinding rate of PLi, :math:`k^\\text{off}_\\text{PLi}`, 
    from fibrin.

    :Units: sec^-1
    :Fortran: kplioff"""

    unbind_rate_tPA_wPLG: Quantity = field(init=False)
    """The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, from fibrin in the presence of PLG.

    Calculated from dissociation constant and binding rate. This is a dependent
    parameter computed in __post_init__().

    :Units: sec^-1
    :Fortran: kaoff12"""

    unbind_rate_tPA_woPLG: Quantity = field(init=False)
    """The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, from fibrin in the absence of PLG.

    Calculated from dissociation constant and binding rate. This is a dependent
    parameter computed in __post_init__().

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
    """The volume fraction of protein in each fiber.

    Calculated from fiber geometry and protofibril packing. This is a dependent
    parameter computed in __post_init__() using equations from Bannish et al. 2017.

    :Units: %
    :Fortran: None"""

    fibrin_conc_per_fiber: Quantity = field(init=False)
    """The concentration of fibrin in each fiber.

    Calculated from fiber geometry and Avogadro's constant. This is a dependent
    parameter computed in __post_init__() using equations from Bannish et al. 2017.

    :Units: micromolar
    :Fortran: None"""

    binding_sites: Quantity = field(init=False)  # int = 427
    """Concentration of binding sites for tPA and plasminogen on each fiber.

    Calculated from the number of doublets in the fiber lattice and the fibrin
    concentration. This is a dependent parameter computed in __post_init__()
    using equations from Bannish et al. 2017.

    :Units: micromolar
    :Fortran: bs"""

    #####################################
    # Model Parameters
    #####################################

    nodes_in_micro_row: int = 7
    """The number of protofibrils in one row of the lattice within a fiber.

    Determines the internal structure of each fiber in the microscale model.
    The total number of protofibrils in a fiber cross-section is nodes_in_micro_row².

    :Units: None
    :Fortran: nodes"""

    snap_proportion: float = 2.0 / 3.0
    """The critical degradation fraction at which a fiber breaks.

    When this proportion of fibrinogen doublets in a fiber have been degraded,
    the fiber is considered to have snapped (completely lysed). Default value
    of 2/3 means the fiber breaks when 67% of doublets are degraded.

    :Units: None
    :Fortran: snap_proportion"""

    #####################################
    # Mechanism Parameters
    #####################################

    micro_simulations: int = 50_000
    """The number of independent trials run in the microscale model.
    
    :Units: None
    :Fortran: simulations"""

    micro_seed: int = 0
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    #####################################
    # Code Parameters
    #####################################

    micro_version: str = "micro_rates"
    """A string identifying which version of the microscale model is being run.

    Used for tracking and logging purposes to distinguish between different
    microscale implementations or parameter sets."""

    micro_log_lvl: int = logging.WARNING
    """The logging level for console output.

    Controls the verbosity of debugging and status information. Uses Python's
    standard logging levels (DEBUG, INFO, WARNING, ERROR, CRITICAL).

    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        """Calculate dependent parameters after initialization.

        This method is automatically called by dataclass.__init__() after all
        independent parameters are set. It computes derived values including:

        - Protofibril radius (2 × fibrinogen radius)
        - Unbinding rates from dissociation constants and binding rates
        - Protein fraction per fiber (from geometric calculations)
        - Fibrin concentration per fiber
        - Binding site concentration

        The method uses object.__setattr__() to set values because the dataclass
        is frozen (immutable).

        Note:
            This method should never be called manually. It runs automatically
            during object construction.
        """

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
                self.nodes_in_micro_row**2
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
                self.nodes_in_micro_row**2
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
            * (self.nodes_in_micro_row - 1)
            / self.nodes_in_micro_row**2
            * self.fibrin_conc_per_fiber,
        )


@dataclass(frozen=True)
class MacroParameters(Parameters):
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

    micro_params: MicroParameters
    """The microscale parameters used to generate input data for this macroscale simulation.

    The macroscale model requires microscale output (binding/unbinding statistics)
    as input. This field holds the microscale parameter set that was used to
    generate that input data."""

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

    # TODO(bpaynter): This value should derive from Microscale Model results
    forced_unbind: float = 0.0852
    """Fraction of times tPA was forced to unbind in microscale model.
    
    :Units: None
    :Fortran: frac_forced"""

    average_bound_time: Quantity = field(init=False)  #  = Q_("27.8 sec")
    """The average time a tPA molecule stays bound to fibrin.

    Calculated as the reciprocal of the unbinding rate in the absence of PLG
    (1 / unbind_rate_tPA_woPLG). This is a dependent parameter computed in
    __post_init__().

    :Units: seconds
    :Fortran: avgwait = 1/kaoff10"""

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
    """The number of rows containing fibrin.

    Calculated as total rows minus empty rows. This is a dependent parameter
    computed in __post_init__().

    :Units: None
    :Fortran: Fhat"""

    empty_rows: int = 29 - 1
    """The number of fibrin-free rows at the top of the grid.

    This represents the depth of the fibrin-free region where tPA molecules
    enter the simulation domain. 
    
    For example, if empty_rows = 28, then rows 0-27 contain no fibers, and
    row 28 is the first row with fibers.

    Note: In the Fortran code, there was no variable for 'The number of empty rows'.
    Instead, the variable `Ffree` gave the 1-indexed location of the first non-empty row.
    In 0-indexing, the number of empty rows is the same as the index of the first 
    non-empty row, but in 1-indexing, this is not the case.
    Thus, `empty_rows` and `Ffree` are not exactly translatable, but
    `empty_rows = Ffree - 1`

    Here is a graphical (rotated) example with four fiber-free rows and 
    five rows of fibrin. 
    In Python, this would give: rows = 9, fiber_rows = 6, and empty_rows = 3.
    In Fortran, the equivalent would be: F = 9, Fhat = 6, and Ffree = 4
    
                      empty_rows
                       ^^^^^^^
    Python indexing:   0  1  2  3  4  5  6  7  8
                       .  .  .  |  |  |  |  |  |
    Fortran indexing:  1  2  3  4  5  6  7  8  9
                                ^
                              Ffree
                                
    :Units: None
    :Fortran: Ffree-1"""

    empty_edges: int = field(init=False)
    """The number of edges without fibrin in the fibrin-free region.

    Also represents the 1-D index of the last edge without fibrin when using
    1-based indexing (for Fortran compatibility). Calculated as full_row × empty_rows.
    This is a dependent parameter computed in __post_init__().

    Note:
        This parameter is primarily for Fortran compatibility and may be
        unnecessary when using 2-D data structures.

    :Units: None
    :Fortran: enoFB"""

    full_row: int = field(init=False)
    """The number of edges in a full row of nodes.

    Calculated as 3 × cols - 1 (right, up, and out edges for each node, except
    the last node which has no right edge). This is a dependent parameter
    computed in __post_init__().

    :Units: None
    :Fortran: None"""

    xz_row: int = field(init=False)
    """The number of horizontal (x) and vertical (z) edges in a row.

    Calculated as 2 × cols - 1 (right and up edges for each node, except the
    last node which has no right edge). This is a dependent parameter computed
    in __post_init__().

    :Units: None
    :Fortran: None"""

    total_edges: int = field(init=False)
    """The total number of edges in the entire grid.

    Calculated as full_row × (rows - 1) + xz_row (full rows for all but the
    last row, which has no up edges). This is a dependent parameter computed
    in __post_init__().

    :Units: None
    :Fortran: num"""

    total_fibers: int = field(init=False)
    """The total number of fibrin fibers in the model.

    Calculated as full_row × (rows - empty_rows - 1) + xz_row, accounting for
    the fibrin-free region and the last row having no up edges. This is a
    dependent parameter computed in __post_init__().

    :Units: None
    :Fortran: None"""

    total_molecules: int = 43074
    """The total number of tPA molecules in the simulation.

    Common values based on physiological concentrations:
        - 43074 corresponds to [tPA] = 0.6 nM
        - 86148 corresponds to [tPA] = 1.2 nM

    :Units: None
    :Fortran: M"""

    moving_probability: float = 0.2
    """The probability that an unbound tPA molecule attempts to move in a timestep.

    This parameter connects the discrete simulation to the continuous diffusion
    equation. It must be small enough to ensure numerical stability and convergence.
    Together with pore size and diffusion coefficient, it determines the timestep
    length via Equation 2.4 (Bannish et al. 2014).

    :Units: None
    :Fortran: q"""

    #####################################
    # Mechanism Parameters
    #####################################

    macro_simulations: int = 10
    """The number of independent simulations to be run
    
    :Units: None
    :Fortran: simulations"""

    total_time: Quantity = Q_("20 min")
    """Total running time for model.
     
    :Units: seconds
    :Fortran: tf"""

    time_step: float = field(init=False)
    """The length of one timestep in the simulation.

    Calculated from the diffusion equation (Equation 2.4, Bannish et al. 2014)
    based on moving probability, pore size, and diffusion coefficient. This is
    a dependent parameter computed in __post_init__().

    :Units: seconds
    :Fortran: tstep"""

    total_time_steps: int = field(init=False)
    """The total number of timesteps in the simulation.

    Calculated as total_time divided by time_step. This is a dependent parameter
    computed in __post_init__().

    :Units: None
    :Fortran: num_t"""

    macro_seed: int = 0  # -2137354075
    """Seed for the random number generator
    
    :Units: None
    :Fortran: seed"""

    state: Tuple[int, int, int, int] = field(init=False)
    """Initial state for the random number generator.

    A 4-tuple of unsigned 32-bit integers where the fourth element is set to
    macro_seed. This is a dependent parameter computed in __post_init__().

    :Units: None
    :Fortran: state"""

    #####################################
    # Data Parameters
    #####################################

    save_interval: Quantity = Q_("10 sec")
    """How often to record data from the model.
    
    :Units: sec
    :Fortran: save_interval"""

    number_of_saves: int = field(init=False)
    """The number of times data will be saved during the simulation.

    Calculated as total_time / save_interval + 1 (one save at the start of each
    interval plus one at the end). This is a dependent parameter computed in
    __post_init__().

    :Units: None
    :Fortran: nplt"""

    #####################################
    # Code Parameters
    #####################################

    macro_version: str = "diffuse_into_and_along"
    """A string identifying which version of the macroscale model is being run.

    Used for tracking and logging purposes to distinguish between different
    macroscale implementations. This string was historically included in data
    filenames by the legacy Fortran code."""

    macro_log_lvl: int = logging.WARNING
    """The logging level for console output.

    Controls the verbosity of debugging and status information. Uses Python's
    standard logging levels (DEBUG, INFO, WARNING, ERROR, CRITICAL).

    :Units: None
    :Fortran: None"""

    duplicate_fortran: bool = False
    """Whether the Python code should replicate Fortran implementation exactly.

    When True, the Python implementation uses the same algorithm sequence and
    random number generation as the legacy Fortran code to produce bit-identical
    results for validation purposes. This mode sacrifices performance for exact
    reproducibility.

    Note:
        This feature requires some modifications to the Fortran code to work.
        See `macro_rng_array.f90`

    :Units: None
    :Fortran: None"""

    processing_library: str = "numpy"
    """Which array processing library the macroscale model should use.

    Options:
        - 'numpy': Standard CPU-based NumPy arrays
        - 'cupy': GPU-accelerated arrays (requires CUDA-compatible GPU)

    The processing library determines whether computations run on CPU or GPU.

    :Units: None
    :Fortran: None"""

    def __post_init__(self):
        """Calculate dependent parameters after initialization.

        This method is automatically called by dataclass.__init__() after all
        independent parameters are set. It computes derived values including:

        - Average bound time (from microscale unbinding rate)
        - Input/output data lists
        - Grid dimension parameters (full_row, xz_row, total_edges, etc.)
        - Fiber count (total_fibers)
        - Empty edge count
        - Time step (from diffusion equation)
        - Total time steps
        - RNG state (initialized with macro_seed)
        - Number of save points

        The method uses object.__setattr__() to set values because the dataclass
        is frozen (immutable).

        Note:
            This method should never be called manually. It runs automatically
            during object construction.
        """
        #
        object.__setattr__(
            self, "average_bound_time", 1.0 / self.micro_params.unbind_rate_tPA_woPLG
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
        object.__setattr__(
            self, "state", (129281, 362436069, 123456789, self.macro_seed)
        )

        # Total saves is one for the start of each 'save_interval' plus one at
        # the end of the run.
        object.__setattr__(
            self, "number_of_saves", int(self.total_time / self.save_interval) + 1
        )
