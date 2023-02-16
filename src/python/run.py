#!/usr/bin/env python3
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

import cProfile
import logging
import os
import sys

from datetime import datetime
from typing import AnyStr

import lysis
import lysis.util


def run(e: lysis.util.Experiment, timestamp: AnyStr):
    if __name__ == "__main__":
        logger = logging.getLogger("lysis")
    else:
        logger = logging.getLogger(__name__)
    logger.info(f"Initialized Experiment '{e.experiment_code}'")
    # p = {
    #     "rows": 12,
    #     "cols": 9,
    #     "empty_rows": 3,
    #     "total_molecules": 430,
    #     "total_time": 10 * 60,
    #     "duplicate_fortran": False,
    # }
    # e.initialize_macro_param(p)
    e.read_file()
    logger.debug(f"With parameters {os.linesep}{e}")
    for file in e.params.output_data:
        filename = lysis.util.default_filenames[file]
        if os.path.isfile(os.path.join(e.os_path, filename)):
            os.remove(os.path.join(e.os_path, filename))
    macro = lysis.MacroscaleRun(e)
    os.makedirs(os.path.join(e.os_path, "macro_pstats"), exist_ok=True)

    filename = "macro_pstats_" + timestamp + ".sts"

    cProfile.runctx(
        "macro.run()",
        globals(),
        locals(),
        filename=os.path.join(e.os_path, "macro_pstats", filename),
    )

    logger.info(f"cProfile stats saved as {filename}.")


def main():
    e = lysis.util.Experiment(r"../../data", experiment_code="2023-01-31-1303")
    timestamp = datetime.now().strftime("%Y-%m-%d-%H%M%S")
    os.makedirs(os.path.join(e.os_path, "log"), exist_ok=True)
    logfile = os.path.join(e.os_path, "log", "lysis-py-" + timestamp + ".log")
    logging.basicConfig(filename=logfile, level=logging.DEBUG)

    formatter = logging.Formatter(
        "%(asctime)s - " "%(name)s - " "%(levelname)s - " "%(message)s"
    )
    logger = logging.getLogger("lysis")
    logger.setLevel(logging.DEBUG)

    stdout = logging.StreamHandler(stream=sys.stdout)
    stdout.setLevel(logging.INFO)
    stdout.setFormatter(formatter)
    logger.addHandler(stdout)

    stderr = logging.StreamHandler(stream=sys.stderr)
    stderr.setLevel(logging.ERROR)
    stderr.setFormatter(formatter)
    logger.addHandler(stderr)

    run(e, timestamp)


if __name__ == "__main__":
    main()
