import inspect
import os
import subprocess

from dataclasses import asdict, dataclass
from typing import AnyStr
from ..experiment import Experiment, Parameters


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

    def exec_command(self):
        arguments = [
            "--expCode",
            self.exp.experiment_code,
            "--inFileCode",
            self.in_file_code,
            "--outFileCode",
            self.out_file_code,
        ]
        params = asdict(self.exp.params)
        sig = inspect.signature(Parameters)
        fortran_names = Parameters.fortran_names()
        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                if fortran_names[key][-2:] == "-1":
                    arguments += ["--" + fortran_names[key][:-2], str(params[key] + 1)]
                else:
                    arguments += ["--" + fortran_names[key], str(params[key])]
        return [self.executable] + arguments

    def run(self):
        command = self.exec_command()
        output_file_name = os.path.join(
            self.exp.os_path,
            os.path.basename(self.executable) + self.out_file_code[:-3] + "txt",
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
