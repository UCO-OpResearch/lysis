"""Utilities for interfacing with Fortran implementations of simulations.

This module provides wrapper classes for executing Fortran-based macroscale
and microscale simulations. It handles parameter conversion, file I/O, and
subprocess management for running compiled Fortran executables.

The main classes:
- FortranMacro: Wrapper for executing the Fortran macroscale simulation
- FortranMicro: Wrapper for executing the Fortran microscale simulation

Both classes convert Python parameters to command-line arguments compatible
with the Fortran programs, manage output files, and handle parallel execution
by splitting random number generator seeds.
"""

import inspect
import os
import subprocess

from dataclasses import asdict, dataclass
from typing import AnyStr

import numpy as np

from pint import Quantity

from .parameters import MacroParameters, MicroParameters
from .run import Run
from .edge_grid import generate_fortran_neighborhood_structure

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


compiler = "ifort"
options = "-r8 -mcmodel medium -traceback"
kiss = "kiss.o"


@dataclass
class FortranMacro:
    """Wrapper for executing the Fortran macroscale simulation.

    This class manages the execution of the compiled Fortran macroscale
    simulation binary. It converts Python Run parameters to command-line
    arguments, generates required input files (neighbor structure), and
    captures output to appropriate files.

    Supports parallel execution by allowing individual simulation instances
    to be run with split random number generator seeds.

    :ivar run: The Run object containing all simulation parameters and data
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran macroscale executable
    :vartype executable: str
    :ivar in_file_code: Code suffix for input files (e.g., for neighbors.dat)
    :vartype in_file_code: str
    :ivar out_file_code: Code suffix for output files (e.g., macro_output.txt)
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs (splits RNG seed if provided)
    :vartype index: int or None
    """

    run: Run = None
    cwd: AnyStr = "."
    #    source: AnyStr = None
    executable: AnyStr = None
    in_file_code: AnyStr = ""
    out_file_code: AnyStr = ""
    index: int = None

    def generate_neighborhoods(self):
        """Generate and write the neighbor structure file for Fortran.

        Pre-calculates the neighborhood structure (8 neighbors for each edge
        in the hexagonal grid) and writes it to a file that the Fortran
        simulation will read. Fortran uses 1-based indexing, so all indices
        are incremented by 1.

        The output file is written to: <run.os_path>/neighbors.dat

        :raises IOError: If unable to write the neighbors.dat file
        """
        fort_neighbors = (
            generate_fortran_neighborhood_structure(
                self.run.macro_params.rows, self.run.macro_params.cols
            )
            + 1
        )
        fort_neighbors.tofile(
            os.path.join(self.run.os_path, "neighbors.dat"), sep=os.linesep
        )

    def exec_command(self):
        """Build the command-line arguments for executing the Fortran binary.

        Converts all macro parameters from the Run object into command-line
        arguments that the Fortran program expects. Handles:

        - Unit conversions for Pint Quantity objects
        - Index adjustments (Fortran uses 1-based indexing)
        - Special scaling factors (e.g., "*100" suffix parameters)
        - RNG seed splitting for parallel runs (if index is set)
        - Micro parameters needed by macro (fiber_radius, binding_sites)

        Only parameters that differ from defaults are included in the
        command line.

        :return: Complete command as a list [executable, arg1, arg2, ...]
            suitable for subprocess.run()
        :rtype: list[str]

        Examples:
            >>> fm = FortranMacro(run=my_run, executable="./macro.exe")
            >>> cmd = fm.exec_command()
            >>> # Returns: ['./macro.exe', '--runCode', 'test01',
            >>>              '--rows', '100', '--cols', '100', ...]
        """
        params = asdict(self.run.macro_params)
        if self.index is not None:
            stream = np.random.SeedSequence(params["macro_seed"])
            seeds = stream.generate_state(params["macro_simulations"])
            params["macro_simulations"] = 1
            params["macro_seed"] = int(np.int32(seeds[self.index]))
            self.out_file_code = self.out_file_code + f"__{self.index:02}"
        arguments = [
            "--runCode",
            self.run.run_code,
            "--inFileCode",
            self.in_file_code,
            "--outFileCode",
            self.out_file_code,
        ]
        sig = inspect.signature(MacroParameters)
        fortran_names = MacroParameters.fortran_names()
        units = MacroParameters.units()
        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                if isinstance(params[key], Quantity):
                    value = params[key].to(units[key]).magnitude
                else:
                    value = params[key]
                if fortran_names[key][-2:] == "-1":
                    arguments += ["--" + fortran_names[key][:-2], str(value + 1)]
                elif fortran_names[key][-4:] == "*100":
                    arguments += [
                        "--" + fortran_names[key][:-4],
                        str(value // 100),
                    ]
                else:
                    arguments += ["--" + fortran_names[key], str(value)]
        arguments += [
            "--radius",
            str(
                self.run.micro_params.fiber_radius.to(
                    MicroParameters.units()["fiber_radius"]
                ).magnitude
            ),
        ]
        arguments += [
            "--bs",
            str(
                self.run.micro_params.binding_sites.to(
                    MicroParameters.units()["binding_sites"]
                ).magnitude
            ),
        ]
        return [self.executable] + arguments

    def exec(self):
        """Execute the Fortran macroscale simulation.

        Performs the complete workflow to run the Fortran simulation:

        1. Generates the neighbor structure file (neighbors.dat)
        2. Builds the command-line arguments
        3. Executes the Fortran binary via subprocess
        4. Redirects stdout to macro<out_file_code>.txt

        The output file is created in the run's output directory.
        Blocks until the Fortran simulation completes.

        :raises FileNotFoundError: If the executable doesn't exist
        :raises subprocess.CalledProcessError: If the Fortran program fails
        :raises IOError: If unable to create output file

        Side effects:
            - Writes neighbors.dat to run.os_path
            - Writes macro<out_file_code>.txt to run.os_path
        """
        self.generate_neighborhoods()
        command = self.exec_command()
        output_file_name = os.path.join(
            self.run.os_path,
            "macro" + self.out_file_code + ".txt",
        )
        with open(output_file_name, "w") as file:
            result = subprocess.run(
                command,
                stdout=file,
                cwd=self.cwd,
            )


@dataclass
class FortranMicro:
    """Wrapper for executing the Fortran microscale simulation.

    This class manages the execution of the compiled Fortran microscale
    simulation binary. It converts Python Run parameters to command-line
    arguments and captures output to appropriate files.

    The microscale simulation generates binding/unbinding statistics that
    are later used by the macroscale simulation. Unlike FortranMacro, this
    doesn't require neighbor structure generation.

    :ivar run: The Run object containing all simulation parameters and data
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran microscale executable
    :vartype executable: str
    :ivar out_file_code: Code suffix for output files (e.g., micro_output.txt)
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs (currently unused/commented out)
    :vartype index: int or None

    Note:
        The index-based seed splitting functionality is currently disabled
        (see commented code in exec_command). All microscale runs use the
        same seed from micro_params.
    """

    run: Run = None
    cwd: AnyStr = "."
    #    source: AnyStr = None
    executable: AnyStr = None
    out_file_code: AnyStr = ""

    def exec_command(self):
        """Build the command-line arguments for executing the Fortran binary.

        Converts all micro parameters from the Run object into command-line
        arguments that the Fortran program expects. Handles:

        - Unit conversions for Pint Quantity objects using m_as()
        - Index adjustments (Fortran uses 1-based indexing, suffix "-1")
        - Special scaling factors (e.g., "*100" suffix parameters)

        Only parameters that differ from defaults are included in the
        command line.

        Note: RNG seed splitting for parallel runs is currently disabled
        (see commented code lines 130-135).

        :return: Complete command as a list [executable, arg1, arg2, ...]
            suitable for subprocess.run()
        :rtype: list[str]

        Examples:
            >>> fm = FortranMicro(run=my_run, executable="./micro.exe")
            >>> cmd = fm.exec_command()
            >>> # Returns: ['./micro.exe', '--runCode', 'test01',
            >>>              '--totalTrials', '1000', ...]
        """
        params = asdict(self.run.micro_params)
        arguments = [
            "--runCode",
            self.run.run_code,
            "--outFileCode",
            self.out_file_code,
        ]
        sig = inspect.signature(MicroParameters)
        fortran_names = MicroParameters.fortran_names()
        units = MicroParameters.units()
        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                if isinstance(params[key], Quantity):
                    params[key] = params[key].m_as(units[key])
                if fortran_names[key][-2:] == "-1":
                    arguments += ["--" + fortran_names[key][:-2], str(params[key] + 1)]
                elif fortran_names[key][-4:] == "*100":
                    arguments += [
                        "--" + fortran_names[key][:-4],
                        str(params[key] // 100),
                    ]
                else:
                    arguments += ["--" + fortran_names[key], str(params[key])]
        return [self.executable] + arguments

    def exec(self):
        """Execute the Fortran microscale simulation.

        Performs the complete workflow to run the Fortran simulation:

        1. Builds the command-line arguments
        2. Executes the Fortran binary via subprocess
        3. Redirects stdout to micro<out_file_code>.txt

        The output file is created in the run's output directory.
        Blocks until the Fortran simulation completes.

        :raises FileNotFoundError: If the executable doesn't exist
        :raises subprocess.CalledProcessError: If the Fortran program fails
        :raises IOError: If unable to create output file

        Side effects:
            - Writes micro<out_file_code>.txt to run.os_path
        """
        command = self.exec_command()
        output_file_name = os.path.join(
            self.run.os_path,
            "micro" + self.out_file_code + ".txt",
        )
        with open(output_file_name, "w") as file:
            result = subprocess.run(
                command,
                stdout=file,
                cwd=self.cwd,
            )

    # def compile(self):
    #     """Compile the Fortran microscale source code (DEPRECATED).
    #
    #     This method attempted to compile the Fortran source using the Intel
    #     Fortran compiler (ifort). It has been disabled because:
    #
    #     - Compilation should be handled separately from execution
    #     - Module loading ("module load oneapi/compiler") doesn't work in subprocess
    #     - Build systems should use proper build tools (Make, CMake, etc.)
    #
    #     Kept as reference for the compilation flags and dependencies.
    #
    #     Original behavior:
    #         - Loaded Intel OneAPI compiler module
    #         - Compiled with: ifort -r8 -mcmodel medium -traceback
    #         - Linked with kiss.o (KISS random number generator)
    #     """
    #     src_path = os.path.dirname(self.source)
    #     bin_path = os.path.dirname(self.executable)
    #     command = [
    #         compiler,
    #         options,
    #         os.path.join(bin_path, kiss),
    #         self.source,
    #         "-o",
    #         self.executable,
    #     ]
    #     subprocess.run("module load oneapi/compiler")
    #     result = subprocess.run(command, stdout=subprocess.PIPE)
    #     print(result.stdout)
