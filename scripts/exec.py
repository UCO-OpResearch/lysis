import cProfile
import logging
import os
import sys

from datetime import datetime
from typing import AnyStr

import lysis
from lysis.config.constants import Q_
from lysis.config.parameters import MacroParameters
from lysis.config.run import Run
from lysis.data.datastore import DataStore


def exec(run: Run, timestamp: AnyStr):
    if __name__ == "__main__":
        logger = logging.getLogger("lysis")
    else:
        logger = logging.getLogger(__name__)
    logger.info(f"Initialized Run '{run.run_code}'")

    # Open DataStore directly from data_root (HDF5 file lives in the root,
    # not inside the run subfolder).
    run.data = DataStore(run.run_code, run.os_data_root)
    logger.info(f"Opened DataStore: {run.data!r}")

    # Load micro_params from the HDF5 file
    run.micro_params = run.data.micro_params

    # Build macro_params from the HDF5 values, overriding total_time to 2 min
    base_macro = run.data.macro_params.to_basedict()
    base_macro["total_time"] = "120 second"
    base_macro["micro_params"] = run.micro_params
    run.macro_params = MacroParameters.parse_from_basedict(base_macro)

    logger.debug(f"With parameters {os.linesep}{run}")
    logger.info(
        f"Simulation: {run.macro_params.total_time_steps:,} timesteps "
        f"({run.macro_params.total_time} total, "
        f"dt={run.macro_params.time_step})"
    )

    macro = lysis.MacroscaleSim(run)
    # os.makedirs(os.path.join(run.os_path, "macro_pstats"), exist_ok=True)

    # filename = "macro_pstats_" + timestamp + ".sts"

    # cProfile.runctx(
    #     "macro.go()",
    #     globals(),
    #     locals(),
    #     filename=os.path.join(run.os_path, "macro_pstats", filename),
    # )
    macro.go()

    # logger.info(f"cProfile stats saved as {filename}.")


def main():
    run = Run(r"data", run_code="2026-02-18-1723")
    timestamp = datetime.now().strftime("%Y-%m-%d-%H%M%S")
    os.makedirs(os.path.join(run.os_path, "log"), exist_ok=True)
    logfile = os.path.join(run.os_path, "log", "lysis-py-" + timestamp + ".log")
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

    exec(run, timestamp)


if __name__ == "__main__":
    main()
