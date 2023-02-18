# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# framework.py
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
"""Framework for experiment components

"""

import json
import os

from abc import ABC, abstractmethod
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import datetime
from typing import Any, List, Mapping, Self, Union

from .parameters import Parameters
from ..util import Const, DataStore, uuid8code, check_current_folder


CONST = Const()


class Experiment:
    """Houses all information about a given experimental run.

    Args:
        data_root: The path of the folder containing datasets

    Attributes:
        experiment_code: The code number of the experiment.
        os_path: The path to the folder containing this experiment's data

    Raises:
        RuntimeError: An invalid data folder was given.
    """

    def __init__(
        self,
        data_root: Union[str, bytes, os.PathLike] = ".",
        experiment_code: str = None,
    ):
        # Check if the data folder path is valid
        if not os.path.isdir(data_root):
            raise RuntimeError("Data folder not found.", data_root)
        self.os_data_root = data_root
        # If no experiment code was given, create a new one from the current
        # date and time.
        if experiment_code is None:
            self.experiment_code = datetime.now().strftime("exp_%Y-%m-%d-%H%M")
        else:
            self.experiment_code = experiment_code

        # Generate the path to the experiment folder and the parameters file
        self.os_path = os.path.join(data_root, str(self.experiment_code))
        os.makedirs(self.os_path, exist_ok=True)

        self.runs = []
        self.tables = []
        self.variations = []

    @staticmethod
    def split_variations(variations):
        if isinstance(variations, list):
            return [
                scenario
                for individual_variation in variations
                for scenario in Experiment.split_variations(individual_variation)
            ]
        elif isinstance(variations, dict):
            scenarios = []
            for key, val in variations.items():
                if isinstance(val, list):
                    for individual_value in val:
                        individual_variation = variations.copy()
                        individual_variation[key] = individual_value
                        scenarios += Experiment.split_variations(individual_variation)
                    break
            if len(scenarios) == 0:
                scenarios = [variations]
            return scenarios


class Scenario:
    def __init__(self, params=None, data=None, short_name="", readme=""):
        self.scenario_id = "scn_" + uuid8code()
        self.short_name = short_name
        self.readme = readme
        # TODO(bpaynter): Check if the parameters are already stored.
        #                 Don't allow parameters to be changed once stored.
        if params is not None:
            self.parameters = Parameters(**params)
        else:
            self.parameters = Parameters()
        if data is not None:
            pass
        else:
            self.data = None

    def __copy__(self):
        pass

    def to_dict(self) -> dict:
        """Returns the internally stored data as a dictionary.

        Does not include system-specific information like paths.
        """
        # Initialize a dictionary of the appropriate parameters
        output = {
            "scenario_id": self.scenario_id,
            "short_name": self.short_name,
            "data_filenames": None,
        }
        # Get the data filenames from the DataStore
        if self.data is not None:
            output["data_filenames"] = self.data.to_dict()
        # Convert the parameters to a dictionary
        if self.parameters is not None:
            output += self.parameters.to_dict()
        return output

    @classmethod
    def load_scenario(cls, path) -> Self:
        """Load the experiment parameters from disk.

        Raises:
            RuntimeError: No parameter file is available for this experiment.
        """
        param_file = os.path.join(path, "scenario.json")
        # Determine whether the parameter file exists for this experiment
        if not os.path.isfile(param_file):
            raise RuntimeError("Scenario parameter file not found.")
        # Open the file
        with open(param_file, "r") as file:
            # Use the JSON library to read in the parameters as a dictionary
            params = json.load(file)
            data_filenames = params.pop("data_filenames", None)
            short_name = params.pop("short_name")
            scenario_id = params.pop("scenario_id")
            if not check_current_folder(path, scenario_id):
                raise RuntimeWarning("Scenario folder incorrectly named.")
        if os.path.isfile(os.path.join(path, "README.rst")):
            with open(os.path.join(path, "README.rst"), "r") as file:
                readme = file.read()
        else:
            readme = ""
        scenario = cls(params, data_filenames, short_name, readme)
        scenario.scenario_id = scenario_id
        return scenario


@dataclass
class Model:
    language = "Python"
    macro_unbind_can_bind = False
    macro_unbind_can_diffuse = False
    micro_unbind_can_bind = False
    micro_unbind_can_diffuse = True
    red_blood_cells = False
    commit = None
    command = None
    seed = None


class Run:
    def __init__(self, exp: Experiment, mod: Model, sc: Scenario, instances: int = 1):
        self.exp = exp
        self.mod = mod
        self.sc = sc
        self.instances = instances
        self.simulations = [Simulation(self, i) for i in range(self.instances)]


class Simulation(ABC):
    def __init__(self, run: Run, instance: int = 0):
        self.run = run
        self.instance = instance

    @abstractmethod
    def run(self):
        pass
