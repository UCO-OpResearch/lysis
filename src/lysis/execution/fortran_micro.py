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
import subprocess

from dataclasses import dataclass
from pathlib import Path

from ..config.parameters import MicroParameters
from ..config.run import Run
from ..dataio.dataspec import dataspec
from ..dataio.datastore import DataStore, HDF5State
from ..dataio.fileops import write_dataset
from .fortran import FortranRunner, MICRO_FORTRAN_DATASPEC_VERSION


#: The eight per-simulation binary datasets in ``microscale_out``.
#:
#: These are the only Fortran-output datasets whose raw ``.dat`` files
#: concatenate meaningfully — every one is shape ``(-1,)`` so
#: :func:`numpy.fromfile` infers length from file size.  The aggregate
#: datasets in ``macroscale_in`` (``tPAleave``, ``tsectPA``, ``lysismat``,
#: ``lenlysisvect``) are histograms / percentile tables produced *per run*
#: and CANNOT be byte-concatenated.
MICROSCALE_BINARY_DATASETS: tuple = (
    "firstPLi",
    "lasttPA",
    "lyscomplete",
    "lysis",
    "PLi",
    "tPA_time",
    "tPAPLiunbd",
    "tPAunbind",
)


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
        # Array path: every sibling task draws from the same fixed-size
        # seed stream so the partition is reproducibility-correct.
        # Legacy single-sim-per-task path (``num_children=None``) keeps
        # the historical ``self.index + 1`` count to preserve bit-for-bit
        # reproducibility of pre-existing single-task artifacts.
        return self.num_children if self.num_children is not None else self.index + 1

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
        num_children: "int | None" = None,
        allow_stale_binary: bool = False,
        skip_binary_verification: bool = False,
        historical_binary_attrs: "dict | None" = None,
        source_stamp: "tuple[str, str] | None" = None,
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
        :param num_children: Total number of sibling array tasks sharing
            this run.  When set together with ``index``,
            :meth:`exec_command` partitions ``micro_simulations`` across
            the children (see :class:`~lysis.execution.fortran.FortranRunner`
            for details).  Defaults to ``None`` (legacy single-sim-per-task).
        :type num_children: int, optional
        :param allow_stale_binary: Downgrade a stamp mismatch between the
            Fortran binary and ``src/fortran/`` from
            :class:`~lysis.tools.provenance.StaleBinaryError` to a
            loud warning, defaults to ``False``.
        :type allow_stale_binary: bool, optional
        :param skip_binary_verification: Bypass the binary↔source
            staleness check entirely (no ``<binary> --version``
            subprocess, no warning banner).  Used by the historical-build
            workflow where the mismatch is intentional and the binary
            may not implement ``--version``.  Defaults to ``False``.
        :type skip_binary_verification: bool, optional
        :param historical_binary_attrs: Pre-computed binary provenance
            dict (from
            :func:`~lysis.tools.provenance.gather_historical_binary_provenance`)
            stamped into the HDF5 file in place of the default
            ``<binary> --version`` query.  Defaults to ``None``.
        :type historical_binary_attrs: dict, optional
        :param source_stamp: Optional pre-computed ``(commit, dirty)``
            tuple for ``src/fortran/``, forwarded to
            :func:`~lysis.tools.provenance.verify_binary_matches_source`.
            Used by Slurm masters that resolve the stamp once before
            copying ``src/`` out to a staging dir, so the child does not
            need to invoke git from outside the checkout.  Defaults to
            ``None`` (in-process git lookup).
        :type source_stamp: tuple[str, str] or None
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

        with DataStore(hdf5_path.stem, str(hdf5_path.parent), mode="r") as ds:
            state = ds.hdf5_state
        if state == HDF5State.MICRO_FILLED:
            raise ValueError(
                f"HDF5 file '{hdf5_path}' already contains microscale "
                "simulation results. The microscale simulation has already "
                f"been run. Expected state: {HDF5State.MICRO_EMPTY.name}."
            )
        if state in (HDF5State.MACRO_EMPTY, HDF5State.MACRO_FILLED, HDF5State.INCONSISTENT):
            raise ValueError(
                f"HDF5 file '{hdf5_path}' already contains macroscale data. "
                "The microscale simulation cannot be run after macroscale "
                "has been initialized. Expected state: "
                f"{HDF5State.MICRO_EMPTY.name}."
            )

        return cls(
            run=run,
            executable=str(Path(executable).resolve()),
            out_file_code=out_file_code,
            index=index,
            num_children=num_children,
            allow_stale_binary=allow_stale_binary,
            skip_binary_verification=skip_binary_verification,
            historical_binary_attrs=historical_binary_attrs,
            source_stamp=source_stamp,
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
    # Slurm-array post-processing
    # ------------------------------------------------------------------

    def concatenate_child_outputs(
        self, data_dir: "Path | str", num_children: int
    ) -> None:
        """Byte-concatenate per-task Fortran outputs into a single set of files.

        The microscale Fortran binary writes its per-simulation results
        as flat 1-D arrays (one record per simulation) to files named
        ``{dataset}{out_file_code}.dat``.  When run as a Slurm array,
        each task ``N`` writes to ``{dataset}{out_file_code}__NN.dat``;
        this method concatenates all ``num_children`` per-task files into
        the unsuffixed filename so :func:`read_data_collection` sees the
        same layout as a single-process run.

        The :file:`params.json` pre-staged by
        :meth:`_write_setup_files` (which carries the *aggregate*
        ``micro_simulations`` count) is left untouched and is the
        authoritative parameter file for the combined collection; any
        per-task ``params.json`` files in *data_dir* are ignored.

        For ``num_children == 1`` no concatenation is needed; the single
        ``__00`` file is renamed to the unsuffixed name and the log file
        is renamed similarly.

        :param data_dir: Directory containing the per-task ``__NN`` output
            files (typically ``{work_dir}/data/{run_code}/``).
        :type data_dir: Path or str
        :param num_children: Number of array tasks that wrote to *data_dir*.
            Must be a positive integer.
        :type num_children: int
        :raises ValueError: If ``num_children`` is less than 1.
        :raises FileNotFoundError: If any expected per-task file is missing.
        """
        if num_children < 1:
            raise ValueError(
                f"num_children must be >= 1, got {num_children}"
            )

        data_dir = Path(data_dir)
        log_prefix = self._log_prefix()

        if num_children == 1:
            # Single-task: just rename rather than concatenate.
            for dataset in MICROSCALE_BINARY_DATASETS:
                src = data_dir / f"{dataset}{self.out_file_code}__00.dat"
                dst = data_dir / f"{dataset}{self.out_file_code}.dat"
                if not src.exists():
                    raise FileNotFoundError(
                        f"Expected per-task file {src} not found "
                        f"(num_children={num_children})"
                    )
                src.rename(dst)
            log_src = data_dir / f"{log_prefix}{self.out_file_code}__00.txt"
            log_dst = data_dir / f"{log_prefix}{self.out_file_code}.txt"
            if log_src.exists():
                log_src.rename(log_dst)
            return

        for dataset in MICROSCALE_BINARY_DATASETS:
            dst = data_dir / f"{dataset}{self.out_file_code}.dat"
            with open(dst, "wb") as out_fh:
                for i in range(num_children):
                    src = data_dir / f"{dataset}{self.out_file_code}__{i:02}.dat"
                    if not src.exists():
                        raise FileNotFoundError(
                            f"Expected per-task file {src} not found "
                            f"(child index {i} of {num_children})"
                        )
                    with open(src, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)

        # read_data_collection also reads the micro log text file; concatenate
        # the per-task logs so the combined dir is a complete collection.
        log_dst = data_dir / f"{log_prefix}{self.out_file_code}.txt"
        with open(log_dst, "wb") as out_fh:
            for i in range(num_children):
                log_src = data_dir / f"{log_prefix}{self.out_file_code}__{i:02}.txt"
                if not log_src.exists():
                    continue
                with open(log_src, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)

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

        # Preflight: confirm the Fortran binary's embedded stamp matches the
        # current src/fortran/ source tree.  Raises StaleBinaryError before
        # we create any output files, unless allow_stale_binary is True (in
        # which case it warns and primes a banner for the log file below).
        self._verify_binary_version()

        data_dir = work_dir / "data" / self.run.run_code
        data_dir.mkdir(parents=True, exist_ok=True)

        # Write params.json for local (non-indexed) runs.
        # When index is set the caller pre-stages the setup files.
        if self.index is None:
            self._write_setup_files(data_dir)

        command = self.exec_command()
        log_file = data_dir / f"micro{self.out_file_code}.txt"
        with open(log_file, "w") as fh:
            self._write_stale_banner(fh)
            subprocess.run(command, stdout=fh, cwd=str(work_dir), check=True)

        return data_dir
