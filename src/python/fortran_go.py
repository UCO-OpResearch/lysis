import argparse
import os

from dataclasses import dataclass

from lysis.util import Run, FortranMacro, FortranMicro

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
    parser.add_argument(
        "executable", type=str, help="The name of the compiled fortran executable."
    )
    parser.add_argument("run_code", type=str)
    parser.add_argument("--in_code", type=str, default=".dat")
    parser.add_argument("--out_code", type=str, default=".dat")
    parser.add_argument(
        "-n", "--index", type=int, help="The index of this simulation in this run."
    )
    parser.add_argument(
        "--cwd", type=str, default=".", help="Root directory of this project."
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    run = Run(os.path.join(args.cwd, "data"), run_code=args.run_code)
    run.read_file()
    if os.path.basename(args.executable)[:5] == "micro":
        fort = FortranMicro(
            run=run,
            cwd=args.cwd,
            executable=args.executable,
            out_file_code=args.out_code,
        )
    elif os.path.basename(args.executable)[:5] == "macro":
        fort = FortranMacro(
            run=run,
            cwd=args.cwd,
            executable=args.executable,
            in_file_code=args.in_code,
            out_file_code=args.out_code,
            index=args.index,
        )
    else:
        raise ValueError(f"Wrong executable type: {args.executable}")
    fort.go()


if __name__ == "__main__":
    main()
