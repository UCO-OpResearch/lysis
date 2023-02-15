import os

from dataclasses import dataclass
from typing import AnyStr

from lysis import Experiment

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass
class FortranMacro:
    exp: Experiment = None
    executable: AnyStr = None
    in_file_code: AnyStr = ".dat"
    out_file_code: AnyStr = ".dat"

    def run(self):
        pass


def main():
    pass


if __name__ == "__main__":
    main()
