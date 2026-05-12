"""Wrapper for executing the Fortran macroscale simulation binary.

Provides :class:`FortranMacro`, which manages execution of the compiled
Fortran macroscale simulation.  It converts Python
:class:`~lysis.config.run.Run` parameters to command-line arguments,
generates the required neighbor-structure input file, and captures binary
stdout to a log file.

The preferred high-level entry point is :meth:`FortranMacro.run_full`.
For step-by-step control use :meth:`FortranMacro.exec_in_workdir` followed
by :meth:`FortranMacro.import_results`.
"""

import subprocess

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import AnyStr

import numpy as np

from ..config.parameters import MacroParameters, MicroParameters
from ..config.run import Run
from ..dataio.dataspec import dataspec
from ..dataio.datastore import DataStore, HDF5State, COMPATIBLE_DATASPEC_VERSION
from ..dataio.fileops import write_dataset, write_data_collection
from ..geometry.edge_grid import generate_fortran_neighborhood_structure
from .fortran import FortranRunner, MACRO_FORTRAN_DATASPEC_VERSION

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass
class FortranMacro(FortranRunner):
    """Wrapper for executing the Fortran macroscale simulation.

    Manages execution of the compiled Fortran macroscale simulation binary.
    Converts Python :class:`~lysis.config.run.Run` parameters to
    command-line arguments, generates the required input files
    (macroscale_in data, ``neighbors.dat``), and captures binary stdout
    to a log file.

    Supports parallel execution by splitting the RNG seed when :attr:`index`
    is set.

    :ivar run: The Run object containing all simulation parameters and data.
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary.
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran macroscale binary.
    :vartype executable: str
    :ivar in_file_code: Code suffix for input files (e.g. ``_PLG2_tPA01_Q2``).
    :vartype in_file_code: str
    :ivar out_file_code: Code suffix for output files (e.g. ``_TF-v__447_798``).
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs.  When set, the RNG seed is
        split so each index produces an independent stream, ``macro_simulations``
        is forced to 1, and ``__{index:02}`` is appended to ``out_file_code``.
    :vartype index: int or None
    """

    in_file_code: AnyStr = ""

    # ------------------------------------------------------------------
    # FortranRunner hooks
    # ------------------------------------------------------------------

    def _get_params(self):
        return self.run.macro_params

    def _get_params_class(self):
        return MacroParameters

    def _seed_field(self) -> str:
        return "macro_seed"

    def _simulations_field(self) -> str:
        return "macro_simulations"

    def _seed_split_count(self, params: dict) -> int:
        return params["macro_simulations"]

    def _base_arguments(self) -> list:
        return [
            "--runCode", self.run.run_code,
            "--inFileCode", self.in_file_code,
            "--outFileCode", self.out_file_code,
        ]

    def _post_arguments(self, params: dict) -> list:
        """Append ``--radius`` and ``--bs`` from micro_params.

        The macroscale binary needs two microscale-level parameters
        (fiber radius and binding sites) that are not part of
        :class:`~lysis.config.parameters.MacroParameters` itself.

        :param params: Unused; present for interface conformance.
        :type params: dict
        :return: CLI arguments for radius and binding sites.
        :rtype: list[str]
        """
        micro_units = MicroParameters.units()
        return [
            "--radius",
            str(self.run.micro_params.fiber_radius.m_as(
                micro_units["fiber_radius"]
            )),
            "--bs",
            str(self.run.micro_params.binding_sites.m_as(
                micro_units["binding_sites"]
            )),
        ]

    def _log_prefix(self) -> str:
        return "macro"

    def _collection_name(self) -> str:
        return "macroscale_out"

    def _fortran_dataspec_version(self) -> str:
        return MACRO_FORTRAN_DATASPEC_VERSION

    def _pre_execute(self) -> None:
        """Generate neighborhoods before the legacy ``execute()`` path runs."""
        self.generate_neighborhoods()

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_hdf5(
        cls,
        hdf5_path: "Path | str",
        executable: str,
        in_file_code: str = "",
        out_file_code: str = "",
        index: "int | None" = None,
        allow_stale_binary: bool = False,
    ) -> "FortranMacro":
        """Construct a :class:`FortranMacro` from an existing HDF5 run file.

        Reads :class:`~lysis.config.parameters.MicroParameters` and
        :class:`~lysis.config.parameters.MacroParameters` from the HDF5 file
        and builds the underlying :class:`~lysis.config.run.Run` automatically.

        :param hdf5_path: Full path to the ``.h5`` file (must contain both
            ``micro_params`` and ``macro_params``).
        :type hdf5_path: Path or str
        :param executable: Path to the compiled Fortran macroscale binary.
        :type executable: str
        :param in_file_code: Input file code suffix, defaults to ``""``.
        :type in_file_code: str, optional
        :param out_file_code: Output file code suffix, defaults to ``""``.
        :type out_file_code: str, optional
        :param index: Parallel run index for seed splitting, defaults to
            ``None``.
        :type index: int, optional
        :param allow_stale_binary: Downgrade a stamp mismatch between the
            Fortran binary and ``src/fortran/`` from
            :class:`~lysis.tools.binary_version.StaleBinaryError` to a
            loud warning, defaults to ``False``.
        :type allow_stale_binary: bool, optional
        :return: Fully configured :class:`FortranMacro` instance.
        :rtype: FortranMacro
        :raises ValueError: If microscale datasets are empty (microscale not yet
            run), if macroscale has not been initialized
            (``initialize_macroscale`` not called), or if macroscale datasets
            are non-empty (macroscale has already been run).
        """
        hdf5_path = Path(hdf5_path)
        run = Run(str(hdf5_path.parent), run_code=hdf5_path.stem)
        run.load_params_from_hdf5()

        with DataStore(hdf5_path.stem, str(hdf5_path.parent), mode="r") as ds:
            state = ds.hdf5_state
        if state in (HDF5State.MICRO_EMPTY, HDF5State.INCONSISTENT):
            raise ValueError(
                f"HDF5 file '{hdf5_path}' does not contain completed "
                "microscale simulation data. The microscale simulation must "
                f"be run before the macroscale simulation. Expected state: "
                f"{HDF5State.MACRO_EMPTY.name}."
            )
        if state == HDF5State.MICRO_FILLED:
            raise ValueError(
                f"HDF5 file '{hdf5_path}' has not been initialized for "
                "macroscale simulation. Call DataStore.initialize_macroscale() "
                f"before running macroscale. Expected state: "
                f"{HDF5State.MACRO_EMPTY.name}."
            )
        if state == HDF5State.MACRO_FILLED:
            raise ValueError(
                f"HDF5 file '{hdf5_path}' already contains macroscale "
                "simulation results. The macroscale simulation has already "
                f"been run. Expected state: {HDF5State.MACRO_EMPTY.name}."
            )

        return cls(
            run=run,
            executable=str(Path(executable).resolve()),
            in_file_code=in_file_code,
            out_file_code=out_file_code,
            index=index,
            allow_stale_binary=allow_stale_binary,
        )

    # ------------------------------------------------------------------
    # Neighborhood generation (legacy path)
    # ------------------------------------------------------------------

    def generate_neighborhoods(self):
        """Generate and write the neighbor-structure file for Fortran.

        Pre-calculates the 8-neighbor structure for every edge in the
        hexagonal grid and writes it to ``neighbors.dat`` in
        ``run.os_path``.  Fortran uses 1-based indexing, so all indices are
        incremented by 1.

        Used by the legacy :meth:`execute` path.  The HDF5-integrated
        :meth:`exec_in_workdir` path generates neighbors via
        :meth:`_generate_macroscale_input_files` instead.

        :raises OSError: If the output file cannot be written.
        """
        fort_neighbors = (
            generate_fortran_neighborhood_structure(
                self.run.macro_params.rows, self.run.macro_params.cols
            )
            + 1  # Convert 0-based (Python) → 1-based (Fortran)
        ).flatten().astype(np.int32)
        write_dataset(
            fort_neighbors,
            self.run.os_path,
            dataspec["v1.99.0"]["macroscale_in"].data["neighbors"],
        )

    # ------------------------------------------------------------------
    # Setup helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _link_setup_files(source_dir: Path, target_dir: Path) -> None:
        """Symlink all non-directory entries from *source_dir* into *target_dir*.

        :param source_dir: Directory containing pre-staged setup files.
        :type source_dir: Path
        :param target_dir: Per-simulation directory to receive the symlinks.
        :type target_dir: Path
        """
        target_dir.mkdir(parents=True, exist_ok=True)
        for entry in source_dir.iterdir():
            if not entry.is_dir():
                link = target_dir / entry.name
                if not link.exists():
                    link.symlink_to(entry)

    @staticmethod
    def _remove_symlinks(target_dir: Path) -> None:
        """Remove all symbolic links from *target_dir*.

        Called after the Fortran binary finishes to leave only its output.

        :param target_dir: Per-simulation directory to clean up.
        :type target_dir: Path
        """
        for entry in target_dir.iterdir():
            if entry.is_symlink():
                entry.unlink()

    def _generate_macroscale_input_files(self, data_dir: Path) -> None:
        """Read microscale output from HDF5 and write v1.99.0 macroscale_in files.

        1. Opens the run's HDF5 DataStore in read-only mode.
        2. Reads the six required microscale datasets as numpy arrays.
        3. Calls :func:`~lysis.dataio.dataconvert.generate_macroscale_in`
           to produce v2.0.0 macroscale_in data (including
           ``edge_grid_neighbors``).
        4. Converts the result to v1.99.0 via
           :func:`~lysis.dataio.dataconvert.convert_data`.
        5. Writes the v1.99.0 text files (``tPAleave``, ``tsectPA``,
           ``lysismat``, ``lenlysisvect``, ``neighbors``) to *data_dir*.

        :param data_dir: Directory to write the macroscale input files into.
        :type data_dir: Path
        """
        # Lazy imports to avoid circular imports
        from ..dataio.dataconvert import generate_macroscale_in, convert_data  # noqa: PLC0415

        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
        dataset_names = [
            "pli_first_time",
            "tpa_leaving_time",
            "fiber_degraded",
            "sim_final_time",
            "tpa_unbound_by_pli",
            "tpa_unbound_kinetic",
        ]

        with DataStore(self.run.run_code, self.run.os_path, mode="r") as ds:
            in_data = {}
            for name in dataset_names:
                ds_spec = micro_spec.data[name]
                in_data[name] = ds._file[ds_spec.data_location][:]

            in_data["params"] = {
                "micro_params": self.run.micro_params.to_basedict(),
                "macro_params": self.run.macro_params.to_basedict(),
            }

        # Generate v2.0.0 macroscale_in (includes edge_grid_neighbors)
        macro_in_v200 = generate_macroscale_in(in_data)

        # Convert to v1.99.0 for Fortran consumption
        macro_in_v199 = convert_data(macro_in_v200, "v2.0.0", "v1.99.0")

        # Write v1.99.0 text files to data_dir
        write_data_collection(
            macro_in_v199,
            str(data_dir),
            [dataspec["v1.99.0"]["macroscale_in"]],
            [self.in_file_code],
        )

    def _write_setup_files(self, data_dir: Path) -> None:
        """Write all shared input files needed before the Fortran binary runs.

        Writes:

        1. ``params.json`` — both ``micro_params`` and ``macro_params``
           (needed by the v1.99.0 import pipeline to resolve dataset shapes).
        2. Macroscale input text files and ``neighbors.dat`` — generated
           from HDF5 microscale output via
           :meth:`_generate_macroscale_input_files`.

        :param data_dir: Directory to write setup files into.
        :type data_dir: Path
        """
        # Write params.json with both micro and macro params
        params_data = {
            "micro_params": self.run.micro_params.to_basedict(),
            "macro_params": self.run.macro_params.to_basedict(),
        }
        write_dataset(
            params_data,
            str(data_dir),
            dataspec[MACRO_FORTRAN_DATASPEC_VERSION]["macroscale_out"].params,
        )

        # Write macroscale_in text files (tPAleave, tsectPA, lysismat,
        # lenlysisvect, neighbors)
        self._generate_macroscale_input_files(data_dir)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def exec_in_workdir(self, work_dir: "Path | str") -> Path:
        """Execute the Fortran macroscale binary for each simulation.

        Each simulation runs in its own isolated subdirectory so that
        concurrent Slurm array tasks never collide.  The Fortran binary is
        invoked with ``--runCode {run_code}/{sim:02}`` so it writes directly to
        ``data/{run_code}/{sim:02}/`` relative to *work_dir*.

        1. Creates ``{work_dir}/data/{run_code}/`` (the shared parent dir).
        2. Writes setup files (``params.json``, macroscale_in text files,
           ``neighbors.dat``) unless :attr:`index` is set (in which case they
           are assumed to be pre-staged by the caller, e.g.
           :func:`~lysis.tools.slurm.submit_macro_slurm_job`).
        3. Determines the simulation range: if :attr:`index` is set, runs only
           that simulation; otherwise runs all ``macro_simulations``.
        4. For each simulation:

           a. Creates ``{data_dir}/{sim:02}/`` and symlinks pre-staged setup
              files from ``{data_dir}/`` into it so Fortran can read its input.
           b. Builds a per-simulation command with
              ``--runCode {run_code}/{sim:02}`` and
              ``--outFileCode {out_file_code}_{sim:02}``.
           c. Runs the Fortran binary with ``cwd=work_dir``, capturing stdout
              to ``{data_dir}/{sim:02}/macro{out_file_code}_{sim:02}.txt``.
           d. Removes the symlinks from the per-sim directory, leaving only
              Fortran output.

        :param work_dir: Working directory for the subprocess.  The binary
            itself must already be present at :attr:`executable` (absolute
            path recommended).
        :type work_dir: Path or str
        :return: Path to ``{work_dir}/data/{run_code}/`` — the directory
            containing all Fortran output subdirectories for this run.
        :rtype: Path
        :raises subprocess.CalledProcessError: If any Fortran invocation exits
            with a non-zero status.
        """
        work_dir = Path(work_dir)

        # Preflight: confirm the Fortran binary's embedded stamp matches the
        # current src/fortran/ source tree.  Raises StaleBinaryError before
        # we create any output files, unless allow_stale_binary is True (in
        # which case it warns and primes a banner for each per-sim log file).
        self._verify_binary_version()

        data_dir = work_dir / "data" / self.run.run_code
        data_dir.mkdir(parents=True, exist_ok=True)

        # Write setup files only for local (non-indexed) runs.
        # When index is set the caller pre-stages the setup files.
        if self.index is None:
            self._write_setup_files(data_dir)

        # Determine simulation range
        n_sims = self.run.macro_params.macro_simulations
        macro_seed = self.run.macro_params.macro_seed
        seeds = np.random.SeedSequence(macro_seed).generate_state(n_sims)

        sims = [self.index] if self.index is not None else list(range(n_sims))

        for sim in sims:
            sim_code = f"{self.out_file_code}_{sim:02}"
            sim_run_code = f"{self.run.run_code}/{sim:02}"

            # Create per-sim subdirectory and symlink setup files into it so
            # Fortran can read its input while writing output to the same dir.
            sim_dir = data_dir / f"{sim:02}"
            self._link_setup_files(data_dir, sim_dir)

            # Build per-simulation params: 1 sim, this sim's seed.
            # ``seeds[sim]`` is a ``np.uint32``; the ``|uint32`` Fortran-tag
            # branch in :meth:`_params_to_arguments` handles the signed-int32
            # cast for the Fortran CLI.
            params = asdict(self.run.macro_params)
            params["macro_simulations"] = 1
            params["macro_seed"] = seeds[sim]

            # Build command with per-sim runCode
            sim_arguments = [
                "--runCode", sim_run_code,
                "--inFileCode", self.in_file_code,
                "--outFileCode", sim_code,
            ]
            sim_arguments += self._params_to_arguments(params)
            sim_arguments += self._post_arguments(params)
            command = [self.executable] + sim_arguments

            # Execute — Fortran writes directly into sim_dir via the runCode path
            log_file = sim_dir / f"macro{sim_code}.txt"
            with open(log_file, "w") as fh:
                self._write_stale_banner(fh)
                subprocess.run(command, stdout=fh, cwd=str(work_dir), check=True)

            # Remove symlinks so only actual Fortran output remains
            self._remove_symlinks(sim_dir)

        return data_dir
