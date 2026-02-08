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

import json
import os
import warnings
from datetime import datetime
from typing import Any, Mapping, Union

from ..data.datastore import DataStore
from ..tools.util import dict_to_formatted_str
from ..config.parameters import MicroParameters, MacroParameters


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2024, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.2"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class Run(object):
    """Houses all information about a given experimental run.

    This object contains:

    * Data location
    * Run parameters (microscale and macroscale)
    * Run data

    It provides methods for:

    * Initializing with default parameters
    * Reading parameters from disk
    * Saving parameters to disk
    * Reading input data from disk
    * Writing result data to disk

    :param data_root: The path of the folder containing datasets
    :type data_root: str, bytes, or os.PathLike
    :param run_code: The code number of the run. This will be the name of the
        folder containing the data specific to this run. This should be a date
        and time in 'YYYY-MM-DD-hhmm' format. If no code is given, one will be
        generated from the current date and time.
    :type run_code: str
    :raises RuntimeError: If an invalid data folder is given.
    """

    def __init__(self, data_root: Union[str, bytes, os.PathLike], run_code: str = None):
        """Initialize a Run object with the given data root and run code.

        :param data_root: The path of the folder containing datasets
        :type data_root: str, bytes, or os.PathLike
        :param run_code: The code number of the run. This will be the name of the
            folder containing the data specific to this run. This should be a date
            and time in 'YYYY-MM-DD-hhmm' format. If no code is given, one will be
            generated from the current date and time.
        :type run_code: str
        :raises RuntimeError: If the data folder is not found.
        """
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
        # self.data = DataStore(self.os_path, default_filenames)

    def __str__(self) -> str:
        """Gives a human-readable, formatted string of the current run's
        parameters.

        :return: A formatted string representation of the run parameters
        :rtype: str
        """
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

        :param params: A dictionary of parameters that differ from the default
            values. For example, {'binding_rate': 10, 'pore_size': 3,}
        :type params: dict, optional
        :raises RuntimeError: If the microscale parameters are not valid.
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

        :param params: A dictionary of parameters that differ from the default
            values. For example, {'binding_rate': 10, 'pore_size': 3,}
        :type params: dict, optional
        :raises RuntimeError: If no microscale parameters are supplied or
            if the macroscale parameters are not valid.
        """
        # The macroscale model is dependent on the parameters and results of the
        # microscale model. If no microscale parameters are supplied, the macroscale
        # model cannot be initialized.
        if self.micro_params is None:
            raise RuntimeError("No Microscale parameters.")
        if params is not None:
            self.macro_params = MacroParameters(self.micro_params, **params)
        else:
            self.macro_params = MacroParameters(micro_params=self.micro_params)

    def to_dict(self) -> dict:
        """Returns the internally stored data as a dictionary.

        Does not include system-specific information like paths.

        :return: A dictionary representation of the run's parameters
        :rtype: dict
        """
        # Initialize a dictionary of the appropriate parameters
        output = {
            "run_code": self.run_code,
            "data_filenames": None,
            "micro_params": None,
            "macro_params": None,
        }
        # Get the data filenames from the DataStore
        # if self.data is not None:
        #     output["data_filenames"] = self.data.to_dict()
        # Convert the Microscale parameters to a dictionary
        if self.micro_params is not None:
            output["micro_params"] = self.micro_params.to_basedict()

        # Convert the Macroscale parameters to a dictionary
        if self.macro_params is not None:
            output["macro_params"] = self.macro_params.to_basedict()
        return output

    def to_file(self) -> None:
        """Stores the run parameters to disk.

        Creates or overwrites the params.json file in the run's data
        folder. This file will contain the current run parameters
        (including any micro- and macroscale parameters) in JSON format.

        :raises RuntimeError: If the parameter file cannot be written or
            if there's an issue with the data serialization.
        """
        with open(self.os_param_file, "w") as file:
            # Convert the internal parameters to a dictionary and then use the
            # JSON module to save to disk.
            json.dump(self.to_dict(), file)

    def read_file(self) -> None:
        """Load the run parameters from disk.

        :raises RuntimeError: If no parameter file is available for this run
            or if the file cannot be read properly.
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
            self.micro_params = MicroParameters.parse_from_basedict(micro_params)
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
            macro_params["micro_params"] = self.micro_params
            self.macro_params = MacroParameters.parse_from_basedict(macro_params)
        else:
            # If there were no parameters in the file, then we leave the
            # object null.
            self.macro_params = None
