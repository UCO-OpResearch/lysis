import argparse
import os

from dataclasses import dataclass
from typing import AnyStr

from lysis.util import Experiment, FortranMacro

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=str, help="The name of the compiled fortran executable.")
    parser.add_argument("exp_code", type=str)
    parser.add_argument("--in_code", type=str, default=".dat")
    parser.add_argument("--out_code", type=str, default=".dat")
    parser.add_argument("-n", "--index", type=int, help="The index of this simulation in this run.")
    parser.add_argument("--cwd", type=str, default=".", help="Root directory of this project.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    e = Experiment(os.path.join(args.cwd, "data"), experiment_code=args.exp_code)
    e.read_file()
    fort = FortranMacro(
        exp=e, 
        cwd=args.cwd, 
        executable=args.executable, 
        in_file_code=args.in_code,
        out_file_code=args.out_code,
        index=args.index,
    )
    fort.run()

if __name__ == "__main__":
    main()
