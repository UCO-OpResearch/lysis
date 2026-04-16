"""Wrapper for executing the Fortran microscale simulation binary.

Provides :class:`FortranMicro`, which manages the complete lifecycle of a
Fortran microscale simulation run:

- Converting Python :class:`~lysis.config.run.Run` parameters to
  command-line arguments.
- Executing the Fortran binary in a temporary working directory.
- Importing the binary output into the HDF5 v2.0.0 format.

The preferred high-level entry point is :meth:`FortranMicro.run_full`.
For step-by-step control use :meth:`FortranMicro.exec_in_workdir` followed
by :meth:`FortranMicro.import_results`.
"""

import subprocess

from dataclasses import dataclass
from pathlib import Path

import h5py

from ..config.parameters import MicroParameters
from ..config.run import Run
from ..dataio.dataspec import dataspec
from ..dataio.datastore import COMPATIBLE_DATASPEC_VERSION
from ..dataio.fileops import write_dataset
from .fortran import FortranRunner, MICRO_FORTRAN_DATASPEC_VERSION

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass
class FortranMicro(FortranRunner):
    """Wrapper for executing the Fortran microscale simulation.

    Manages execution of the compiled Fortran microscale simulation binary.
    Converts Python :class:`~lysis.config.run.Run` parameters to
    command-line arguments, captures output to appropriate files, and (via
    :meth:`import_results`) converts the Fortran output into the project's
    HDF5 v2.0.0 format.

    The Fortran binary writes its data files to ``data/{run_code}/`` relative
    to the working directory it is given (``cwd``).  All high-level methods
    here account for this path structure.

    :ivar run: The Run object containing all simulation parameters and data.
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary.
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran microscale binary.
    :vartype executable: str
    :ivar out_file_code: Code suffix for output files (e.g. ``_PLG2_tPA01_Q4``).
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs.  When set, the RNG seed is
        split so each index produces an independent stream, ``micro_simulations``
        is forced to 1, and ``__{index:02}`` is appended to ``out_file_code``.
    :vartype index: int or None
    """

    # ------------------------------------------------------------------
    # FortranRunner hooks
    # ------------------------------------------------------------------

    def _get_params(self):
        return self.run.micro_params

    def _get_params_class(self):
        return MicroParameters

    def _seed_field(self) -> str:
        return "micro_seed"

    def _simulations_field(self) -> str:
        return "micro_simulations"

    def _seed_split_count(self, params: dict) -> int:
        return self.index + 1

    def _base_arguments(self) -> list:
        return [
            "--runCode", self.run.run_code,
            "--outFileCode", self.out_file_code,
        ]

    def _log_prefix(self) -> str:
        return "micro"

    def _collection_name(self) -> str:
        return "microscale_out"

    def _fortran_dataspec_version(self) -> str:
        return MICRO_FORTRAN_DATASPEC_VERSION

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_hdf5(
        cls,
        hdf5_path: "Path | str",
        executable: str,
        out_file_code: str = "",
        index: "int | None" = None,
    ) -> "FortranMicro":
        """Construct a :class:`FortranMicro` from an existing HDF5 run file.

        Reads :class:`~lysis.config.parameters.MicroParameters` from the
        HDF5 file and builds the underlying :class:`~lysis.config.run.Run`
        automatically.

        :param hdf5_path: Full path to the ``.h5`` file (must contain
            ``micro_params``).
        :type hdf5_path: Path or str
        :param executable: Path to the compiled Fortran microscale binary.
        :type executable: str
        :param out_file_code: Output file code suffix, defaults to ``""``.
        :type out_file_code: str, optional
        :param index: Parallel run index for seed splitting, defaults to
            ``None``.
        :type index: int, optional
        :return: Fully configured :class:`FortranMicro` instance.
        :rtype: FortranMicro
        :raises RuntimeError: If the HDF5 file's directory is not found.
        :raises ValueError: If the HDF5 file has an incompatible dataspec, if
            macroscale data is already present, or if microscale datasets are
            non-empty (simulation has already been run).
        """
        hdf5_path = Path(hdf5_path)
        run = Run(str(hdf5_path.parent), run_code=hdf5_path.stem)
        run.load_params_from_hdf5()

        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
        with h5py.File(hdf5_path, "r") as f:
            if "macro_data" in f:
                raise ValueError(
                    f"HDF5 file '{hdf5_path}' already contains macroscale data. "
                    "The microscale simulation cannot be run after macroscale "
                    "has been initialized. Expected file state: empty microscale "
                    "datasets, no macroscale data (post init-experiment)."
                )
            check_loc = micro_spec.data["tpa_leaving_time"].data_location
            if check_loc in f and f[check_loc].shape[0] != 0:
                raise ValueError(
                    f"HDF5 file '{hdf5_path}' already contains microscale "
                    "simulation results. The microscale simulation has already "
                    "been run. Expected file state: empty microscale datasets, "
                    "no macroscale data (post init-experiment)."
                )

        return cls(
            run=run,
            executable=str(Path(executable).resolve()),
            out_file_code=out_file_code,
            index=index,
        )

    # ------------------------------------------------------------------
    # Setup helpers
    # ------------------------------------------------------------------

    def _write_setup_files(self, data_dir: Path) -> None:
        """Write ``params.json`` needed by the v1.99.0 import pipeline.

        :param data_dir: Directory to write the file into.
        :type data_dir: Path
        """
        params_data = {"micro_params": self.run.micro_params.to_basedict()}
        write_dataset(
            params_data,
            str(data_dir),
            dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"].params,
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def exec_in_workdir(self, work_dir: "Path | str") -> Path:
        """Execute the Fortran binary using *work_dir* as the working directory.

        The Fortran code always writes its binary output to
        ``data/{run_code}/`` relative to its working directory.  This method:

        1. Creates ``{work_dir}/data/{run_code}/``.
        2. Writes ``params.json`` there (needed by the v1.99.0 import pipeline
           to resolve shapes) unless :attr:`index` is set, in which case the
           caller is assumed to have pre-staged it.
        3. Executes the binary with ``cwd=work_dir``, redirecting stdout
           (the log file) to ``{work_dir}/data/{run_code}/micro{out_code}.txt``
           so that all import inputs are co-located.

        :param work_dir: Working directory for the subprocess.  Must be on a
            filesystem with enough space for the binary output.  The binary
            itself must already be present at :attr:`executable` (absolute
            path recommended).
        :type work_dir: Path or str
        :return: Path to ``{work_dir}/data/{run_code}/`` — the directory
            containing all Fortran output for this run.
        :rtype: Path
        :raises subprocess.CalledProcessError: If the Fortran program exits
            with a non-zero status.
        """
        work_dir = Path(work_dir)
        data_dir = work_dir / "data" / self.run.run_code
        data_dir.mkdir(parents=True, exist_ok=True)

        # Write params.json for local (non-indexed) runs.
        # When index is set the caller pre-stages the setup files.
        if self.index is None:
            self._write_setup_files(data_dir)

        command = self.exec_command()
        log_file = data_dir / f"micro{self.out_file_code}.txt"
        with open(log_file, "w") as fh:
            subprocess.run(command, stdout=fh, cwd=str(work_dir), check=True)

        return data_dir
