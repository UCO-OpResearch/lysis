import inspect
import os
import subprocess

from dataclasses import asdict, dataclass
from typing import AnyStr

import numpy as np

from pint import Quantity

from .parameters import Experiment, MacroParameters, MicroParameters
from .edge_grid import EdgeGrid

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


compiler = "ifort"
options = "-r8 -mcmodel medium -traceback"
kiss = "kiss.o"


@dataclass
class FortranMacro:
    exp: Experiment = None
    cwd: AnyStr = "."
    #    source: AnyStr = None
    executable: AnyStr = None
    in_file_code: AnyStr = ".dat"
    out_file_code: AnyStr = ".dat"
    index: int = None
    
    def generate_neighborhoods(self):
        fort_neighbors = EdgeGrid.generate_fortran_neighborhood_structure(self.exp) + 1
        fort_neighbors.tofile(os.path.join(self.exp.os_path, "neighbors.dat"), sep=os.linesep)

    def exec_command(self):
        params = asdict(self.exp.macro_params)
        if self.index is not None:
            stream = np.random.SeedSequence(params['seed'])
            seeds = stream.generate_state(params['total_trials'])
            params['total_trials'] = 1
            params['seed'] = int(np.int32(seeds[self.index]))
            self.out_file_code = self.out_file_code[:-4] + f"_{self.index:02}" + self.out_file_code[-4:]
        arguments = [
            "--runCode",
            self.exp.experiment_code,
            "--inFileCode",
            self.in_file_code,
            "--outFileCode",
            self.out_file_code,
        ]
        sig = inspect.signature(MacroParameters)
        fortran_names = MacroParameters.fortran_names()
        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                if fortran_names[key][-2:] == "-1":
                    arguments += ["--" + fortran_names[key][:-2], str(params[key] + 1)]
                elif fortran_names[key][-4:] == "*100":
                    arguments += ["--" + fortran_names[key][:-4], str(params[key] // 100)]
                else:
                    arguments += ["--" + fortran_names[key], str(params[key])]
        return [self.executable] + arguments

    def run(self):
        self.generate_neighborhoods()
        command = self.exec_command()
        output_file_name = os.path.join(
            self.exp.os_path,
            "macro" + self.out_file_code[:-3] + "txt",
        )
        with open(output_file_name, "w") as file:
            result = subprocess.run(
                command,
                stdout=file,
                cwd=self.cwd,
            )


@dataclass
class FortranMicro:
    exp: Experiment = None
    cwd: AnyStr = "."
    #    source: AnyStr = None
    executable: AnyStr = None
    out_file_code: AnyStr = ".dat"
    index: int = None
    
    def exec_command(self):
        params = asdict(self.exp.micro_params)
        # if self.index is not None:
        #     stream = np.random.SeedSequence(params['seed'])
        #     seeds = stream.generate_state(params['total_trials'])
        #     params['total_trials'] = 1
        #     params['seed'] = int(np.int32(seeds[self.index]))
        #     self.out_file_code = self.out_file_code[:-4] + f"_{self.index:02}" + self.out_file_code[-4:]
        arguments = [
            "--runCode",
            self.exp.experiment_code,
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
                    arguments += ["--" + fortran_names[key][:-4], str(params[key] // 100)]
                else:
                    arguments += ["--" + fortran_names[key], str(params[key])]
        return [self.executable] + arguments

    def run(self):
        command = self.exec_command()
        output_file_name = os.path.join(
            self.exp.os_path,
            "micro" + self.out_file_code[:-3] + "txt",
        )
        with open(output_file_name, "w") as file:
            result = subprocess.run(
                command,
                stdout=file,
                cwd=self.cwd,
            )

    # def compile(self):
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
