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

import shutil
import tempfile

from dataclasses import dataclass
from pathlib import Path

from ..config.parameters import MicroParameters
from ..config.run import Run
from ..dataio.dataspec import dataspec
from ..dataio.datastore import DataStore
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
        :raises ValueError: If the HDF5 file has an incompatible dataspec.
        """
        hdf5_path = Path(hdf5_path)
        run = Run(str(hdf5_path.parent), run_code=hdf5_path.stem)
        run.load_params_from_hdf5()
        return cls(
            run=run,
            executable=str(Path(executable).resolve()),
            out_file_code=out_file_code,
            index=index,
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def execute(self) -> None:
        """Execute the Fortran microscale simulation (legacy path).

        Runs the binary with ``cwd=self.cwd`` and writes stdout to
        ``micro{out_file_code}.txt`` inside ``run.os_path``.

        For HDF5-integrated workflows prefer :meth:`exec_in_workdir` /
        :meth:`run_full` instead.

        :raises subprocess.CalledProcessError: If the Fortran program exits
            with a non-zero status.
        """
        super().execute()

    def exec_in_workdir(self, work_dir: "Path | str") -> Path:
        """Execute the Fortran binary using *work_dir* as the working directory.

        The Fortran code always writes its binary output to
        ``data/{run_code}/`` relative to its working directory.  This method:

        1. Creates ``{work_dir}/data/{run_code}/``.
        2. Writes a ``params.json`` there (needed by the v1.99.0 import
           pipeline to resolve shapes when reading binary files back).
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
        import subprocess

        work_dir = Path(work_dir)
        data_dir = work_dir / "data" / self.run.run_code
        data_dir.mkdir(parents=True, exist_ok=True)

        # Write params.json so the v1.99.0 read pipeline can resolve parameters.
        params_data = {"micro_params": self.run.micro_params.to_basedict()}
        write_dataset(
            params_data,
            str(data_dir),
            dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"].params,
        )

        command = self.exec_command()
        log_file = data_dir / f"micro{self.out_file_code}.txt"
        with open(log_file, "w") as fh:
            subprocess.run(command, stdout=fh, cwd=str(work_dir), check=True)

        return data_dir

    # ------------------------------------------------------------------
    # Import
    # ------------------------------------------------------------------

    def import_results(
        self,
        data_dir: "Path | str",
        *,
        keep_on_failure: bool = False,
        keep_tmpdir: bool = False,
    ) -> None:
        """Convert Fortran output in *data_dir* and write it to the run's HDF5 file.

        Opens the run's HDF5 file (``{run.os_path}/{run.run_code}.h5``) in
        append mode and delegates to
        :meth:`~lysis.dataio.datastore.DataStore.import_collection` to read
        the Fortran microscale output (version
        :data:`~lysis.execution.fortran.MICRO_FORTRAN_DATASPEC_VERSION`),
        convert it to HDF5 v2.0.0, and populate the empty ``microscale_out``
        datasets.

        The target HDF5 file should already contain ``micro_params`` and
        empty ``microscale_out`` datasets (i.e.
        :meth:`~lysis.dataio.datastore.DataStore.create` should have been
        called beforehand).

        Cleanup behaviour:

        - On **success**: *data_dir* is removed unless ``keep_tmpdir=True``.
        - On **failure**: *data_dir* is preserved when
          ``keep_on_failure=True`` *or* ``keep_tmpdir=True``; otherwise it
          is removed.

        :param data_dir: Directory containing the Fortran binary output files
            (``lysis{code}.dat``, ``micro{code}.txt``, ``params.json``, …).
            Typically ``{work_dir}/data/{run_code}/``.
        :type data_dir: Path or str
        :param keep_on_failure: Preserve *data_dir* if import raises an
            exception, defaults to ``False``.
        :type keep_on_failure: bool, optional
        :param keep_tmpdir: Always preserve *data_dir*, defaults to ``False``.
        :type keep_tmpdir: bool, optional
        :raises Exception: Re-raises any exception from the import pipeline
            after applying the cleanup policy.
        """
        data_dir = Path(data_dir)

        try:
            with DataStore(self.run.run_code, self.run.os_path, mode="a") as ds:
                ds.import_collection(
                    "microscale_out",
                    MICRO_FORTRAN_DATASPEC_VERSION,
                    str(data_dir),
                    [self.out_file_code],
                )

            if not keep_tmpdir:
                shutil.rmtree(data_dir)

        except Exception:
            if not keep_on_failure and not keep_tmpdir:
                shutil.rmtree(data_dir, ignore_errors=True)
            raise

    # ------------------------------------------------------------------
    # Full workflow
    # ------------------------------------------------------------------

    def run_full(
        self,
        *,
        keep_tmpdir: bool = False,
    ) -> None:
        """Execute the complete HDF5-integrated workflow.

        Creates a uniquely named temporary working directory (via
        :func:`tempfile.mkdtemp` with an explicit *dir* argument so that
        ``/tmp`` — which is a ramdisk on many HPC nodes — is never used),
        runs the Fortran binary, imports the results into the run's HDF5 file,
        and cleans up.

        The temporary directory is placed in :attr:`run.os_path
        <lysis.config.run.Run.os_path>` (the run's data directory) to avoid
        using ramdisk storage.

        Cleanup policy:

        - **Success**: temporary directory removed unless ``keep_tmpdir=True``.
        - **Failure**: temporary directory is preserved for debugging
          (``import_results`` is called with ``keep_on_failure=True``); the
          outer ``except`` block removes it only when ``keep_tmpdir=False``.

        :param keep_tmpdir: Always preserve the temporary directory,
            defaults to ``False``.
        :type keep_tmpdir: bool, optional
        :raises subprocess.CalledProcessError: If the Fortran binary fails.
        :raises Exception: Re-raises any exception from the import pipeline.
        """
        tmpdir = Path(
            tempfile.mkdtemp(
                prefix=f"lysis-micro-{self.run.run_code}-",
                dir=self.run.os_path,
            )
        )
        try:
            data_dir = self.exec_in_workdir(tmpdir)
            self.import_results(
                data_dir,
                keep_on_failure=True,
                keep_tmpdir=keep_tmpdir,
            )
        except Exception:
            if not keep_tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)
            raise
        else:
            if not keep_tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)
