"""Utilities for interfacing with Fortran implementations of simulations.

This module provides wrapper classes for executing Fortran-based macroscale
and microscale simulations. It handles parameter conversion, file I/O, and
subprocess management for running compiled Fortran executables.

The main classes:
- FortranMacro: Wrapper for executing the Fortran macroscale simulation
- FortranMicro: Wrapper for executing the Fortran microscale simulation

Both classes convert Python parameters to command-line arguments compatible
with the Fortran programs, manage output files, and handle parallel execution
by splitting random number generator seeds.
"""

import inspect
import json
import os
import shutil
import subprocess
import tempfile

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import AnyStr

import h5py
import numpy as np

from pint import Quantity

from ..config.parameters import MacroParameters, MicroParameters
from ..config.run import Run
from ..geometry.edge_grid import generate_fortran_neighborhood_structure

#: Dataspec version produced by the Fortran microscale binary.
#: Update this constant when the binary's output format changes.
MICRO_FORTRAN_DATASPEC_VERSION: str = "v1.99.0"

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# Intel Fortran compiler executable
compiler = "ifort"
# Compiler flags: -r8 (promote reals to 8 bytes), -mcmodel medium (large memory model),
# -traceback (enable runtime stack trace on errors)
options = "-r8 -mcmodel medium -traceback"
# Pre-compiled KISS random number generator object file (must be linked with Fortran code)
kiss = "kiss.o"


@dataclass
class FortranMacro:
    """Wrapper for executing the Fortran macroscale simulation.

    This class manages the execution of the compiled Fortran macroscale
    simulation binary. It converts Python Run parameters to command-line
    arguments, generates required input files (neighbor structure), and
    captures output to appropriate files.

    Supports parallel execution by allowing individual simulation instances
    to be run with split random number generator seeds.

    :ivar run: The Run object containing all simulation parameters and data
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran macroscale executable
    :vartype executable: str
    :ivar in_file_code: Code suffix for input files (e.g., for neighbors.dat)
    :vartype in_file_code: str
    :ivar out_file_code: Code suffix for output files (e.g., macro_output.txt)
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs (splits RNG seed if provided)
    :vartype index: int or None
    """

    run: Run = None
    cwd: AnyStr = "."  # Current working directory, defaults to current location
    #    source: AnyStr = None  # Commented out: source path no longer needed (compile separately)
    executable: AnyStr = None  # Path to compiled Fortran executable
    in_file_code: AnyStr = ""  # Suffix for input files (e.g., "_PLG2_tPA01_Q2")
    out_file_code: AnyStr = ""  # Suffix for output files (e.g., "_TF-v__447_798")
    index: int = None  # Parallel run index (if set, splits RNG seed)

    def generate_neighborhoods(self):
        """Generate and write the neighbor structure file for Fortran.

        Pre-calculates the neighborhood structure (8 neighbors for each edge
        in the hexagonal grid) and writes it to a file that the Fortran
        simulation will read. Fortran uses 1-based indexing, so all indices
        are incremented by 1.

        The output file is written to: <run.os_path>/neighbors.dat

        :raises IOError: If unable to write the neighbors.dat file
        """
        # Generate neighborhood structure using Fortran indexing conventions
        fort_neighbors = (
            generate_fortran_neighborhood_structure(
                self.run.macro_params.rows, self.run.macro_params.cols
            )
            + 1  # Add 1 to convert from Python's 0-based to Fortran's 1-based indexing
        )
        # Write the neighbor array to file, one neighbor index per line
        fort_neighbors.tofile(
            os.path.join(self.run.os_path, "neighbors.dat"), sep=os.linesep
        )

    def exec_command(self):
        """Build the command-line arguments for executing the Fortran binary.

        Converts all macro parameters from the Run object into command-line
        arguments that the Fortran program expects. Handles:

        - Unit conversions for Pint Quantity objects
        - Index adjustments (Fortran uses 1-based indexing)
        - Special scaling factors (e.g., "*100" suffix parameters)
        - RNG seed splitting for parallel runs (if index is set)
        - Micro parameters needed by macro (fiber_radius, binding_sites)

        Only parameters that differ from defaults are included in the
        command line.

        :return: Complete command as a list [executable, arg1, arg2, ...]
            suitable for subprocess.run()
        :rtype: list[str]

        Examples:
            >>> fm = FortranMacro(run=my_run, executable="./macro.exe")
            >>> cmd = fm.exec_command()
            >>> # Returns: ['./macro.exe', '--runCode', 'test01',
            >>>              '--rows', '100', '--cols', '100', ...]
        """
        # Convert macro parameters dataclass to dictionary for easier manipulation
        params = asdict(self.run.macro_params)
        # If running in parallel mode, split the RNG seed
        if self.index is not None:
            # Create a seed sequence from the base seed
            stream = np.random.SeedSequence(params["macro_seed"])
            # Generate independent seeds for each parallel run
            seeds = stream.generate_state(params["macro_simulations"])
            # This instance will run only one simulation
            params["macro_simulations"] = 1
            # Use the seed specific to this index, cast to int32 to match Fortran INTEGER*4
            params["macro_seed"] = int(np.int32(seeds[self.index]))
            # Append index to output file code for unique output files
            self.out_file_code = self.out_file_code + f"__{self.index:02}"

        # Start building the command-line argument list
        arguments = [
            "--runCode",
            self.run.run_code,  # Unique identifier for this run
            "--inFileCode",
            self.in_file_code,  # Input file suffix
            "--outFileCode",
            self.out_file_code,  # Output file suffix
        ]

        # Get parameter metadata for conversions
        sig = inspect.signature(MacroParameters)  # Default values
        fortran_names = MacroParameters.fortran_names()  # Map Python names to Fortran names
        units = MacroParameters.units()  # Expected units for each parameter

        # Convert each parameter to command-line argument
        for key in sig.parameters:
            # Only include parameters that have Fortran equivalents and differ from defaults
            if key in fortran_names and params[key] != sig.parameters[key].default:
                # Handle unit conversions for Pint Quantity objects
                if isinstance(params[key], Quantity):
                    # Convert to expected units and extract magnitude (numeric value)
                    value = params[key].to(units[key]).magnitude
                else:
                    value = params[key]

                # Handle special parameter name suffixes that indicate needed transformations
                if fortran_names[key][-2:] == "-1":
                    # Suffix "-1" means: convert 0-based index to 1-based (Python → Fortran)
                    # Strip the suffix and add 1 to the value
                    arguments += ["--" + fortran_names[key][:-2], str(value + 1)]
                elif fortran_names[key][-4:] == "*100":
                    # Suffix "*100" means: parameter is stored as x100 internally
                    # Strip the suffix and divide value by 100
                    arguments += [
                        "--" + fortran_names[key][:-4],
                        str(value // 100),
                    ]
                else:
                    # No special handling needed, pass value as-is
                    arguments += ["--" + fortran_names[key], str(value)]

        # Add micro parameters needed by macro (fiber geometry and binding sites)
        # These come from micro_params but are needed by the macro simulation
        arguments += [
            "--radius",
            str(
                self.run.micro_params.fiber_radius.to(
                    MicroParameters.units()["fiber_radius"]
                ).magnitude  # Convert to expected units
            ),
        ]
        arguments += [
            "--bs",  # bs = binding sites
            str(
                self.run.micro_params.binding_sites.to(
                    MicroParameters.units()["binding_sites"]
                ).magnitude  # Convert to expected units
            ),
        ]
        # Return complete command: [executable, arg1, val1, arg2, val2, ...]
        return [self.executable] + arguments

    def exec(self):
        """Execute the Fortran macroscale simulation.

        Performs the complete workflow to run the Fortran simulation:

        1. Generates the neighbor structure file (neighbors.dat)
        2. Builds the command-line arguments
        3. Executes the Fortran binary via subprocess
        4. Redirects stdout to macro<out_file_code>.txt

        The output file is created in the run's output directory.
        Blocks until the Fortran simulation completes.

        :raises FileNotFoundError: If the executable doesn't exist
        :raises subprocess.CalledProcessError: If the Fortran program fails
        :raises IOError: If unable to create output file

        Side effects:
            - Writes neighbors.dat to run.os_path
            - Writes macro<out_file_code>.txt to run.os_path
        """
        # Generate the neighbor structure file that Fortran will read
        self.generate_neighborhoods()
        # Build the complete command-line arguments
        command = self.exec_command()
        # Construct output file path
        output_file_name = os.path.join(
            self.run.os_path,
            "macro" + self.out_file_code + ".txt",
        )
        # Execute Fortran binary and redirect stdout to file
        with open(output_file_name, "w") as file:
            result = subprocess.run(
                command,
                stdout=file,  # Capture all Fortran output to file
                cwd=self.cwd,  # Execute from specified working directory
            )


@dataclass
class FortranMicro:
    """Wrapper for executing the Fortran microscale simulation.

    This class manages the execution of the compiled Fortran microscale
    simulation binary. It converts Python :class:`~lysis.config.run.Run`
    parameters to command-line arguments, captures output to appropriate
    files, and (via :meth:`import_results`) converts the Fortran output into
    the project's HDF5 v2.0.0 format.

    The Fortran binary writes its data files to ``data/{run_code}/`` relative
    to the working directory it is given (``cwd``).  All high-level methods
    here account for this path structure.

    :ivar run: The Run object containing all simulation parameters and data
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran microscale executable
    :vartype executable: str
    :ivar out_file_code: Code suffix for output files (e.g., ``_PLG2_tPA01_Q4``)
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs.  When set, the RNG seed is
        split so each index produces an independent stream and
        ``micro_simulations`` is forced to 1.  Append ``__{index:02}`` to
        ``out_file_code`` automatically.
    :vartype index: int or None
    """

    run: Run = None
    cwd: AnyStr = "."
    executable: AnyStr = None
    out_file_code: AnyStr = ""
    index: int = None

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
        return cls(run=run, executable=str(executable),
                   out_file_code=out_file_code, index=index)

    # ------------------------------------------------------------------
    # Command building
    # ------------------------------------------------------------------

    def exec_command(self) -> list:
        """Build the command-line arguments for executing the Fortran binary.

        Converts all micro parameters from the Run object into command-line
        arguments that the Fortran program expects. Handles:

        - RNG seed splitting when ``self.index`` is set (mirrors
          :class:`FortranMacro` behaviour)
        - Unit conversions for Pint Quantity objects using ``m_as()``
        - Index adjustments (Fortran uses 1-based indexing, suffix ``"-1"``)
        - Special scaling factors (suffix ``"*100"``)

        Only parameters that differ from their defaults are included.

        .. note::
            When ``self.index`` is not ``None`` this method mutates
            ``self.out_file_code`` by appending ``__{index:02}``.  Calling it
            more than once will keep appending — call it exactly once per
            execution.

        :return: Complete command as ``[executable, arg1, val1, ...]``
            suitable for :func:`subprocess.run`.
        :rtype: list[str]
        """
        params = asdict(self.run.micro_params)

        # Seed splitting for parallel (multi-node) runs
        if self.index is not None:
            stream = np.random.SeedSequence(params["micro_seed"])
            seeds = stream.generate_state(self.index + 1)
            params["micro_simulations"] = 1
            params["micro_seed"] = int(np.int32(seeds[self.index]))
            self.out_file_code = self.out_file_code + f"__{self.index:02}"

        arguments = [
            "--runCode",
            self.run.run_code,
            "--outFileCode",
            self.out_file_code,
        ]

        sig = inspect.signature(MicroParameters)
        fortran_names = MicroParameters.fortran_names()
        units = MicroParameters.units()

        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                if isinstance(params[key], Quantity):
                    params[key] = params[key].m_as(units[key])

                if fortran_names[key][-2:] == "-1":
                    arguments += ["--" + fortran_names[key][:-2], str(params[key] + 1)]
                elif fortran_names[key][-4:] == "*100":
                    arguments += [
                        "--" + fortran_names[key][:-4],
                        str(params[key] // 100),
                    ]
                else:
                    arguments += ["--" + fortran_names[key], str(params[key])]

        return [self.executable] + arguments

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def exec(self) -> None:
        """Execute the Fortran microscale simulation (legacy path).

        Runs the binary with ``cwd=self.cwd`` and writes stdout to
        ``micro{out_file_code}.txt`` inside ``run.os_path``.

        For HDF5-integrated workflows prefer :meth:`exec_in_workdir` /
        :meth:`run_full` instead.

        :raises subprocess.CalledProcessError: If the Fortran program exits
            with a non-zero status.
        """
        command = self.exec_command()
        output_file_name = os.path.join(
            self.run.os_path,
            "micro" + self.out_file_code + ".txt",
        )
        with open(output_file_name, "w") as file:
            subprocess.run(command, stdout=file, cwd=self.cwd)

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
            itself must already be present at ``self.executable`` (absolute
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

        # Write params.json so the v1.99.0 read pipeline can resolve parameters.
        params_data = {"micro_params": self.run.micro_params.to_basedict()}
        with open(data_dir / "params.json", "w") as fh:
            json.dump(params_data, fh, indent=4, default=str)

        command = self.exec_command()
        log_file = data_dir / f"micro{self.out_file_code}.txt"
        with open(log_file, "w") as fh:
            subprocess.run(command, stdout=fh, cwd=str(work_dir), check=True)

        return data_dir

    # ------------------------------------------------------------------
    # Import
    # ------------------------------------------------------------------

    @staticmethod
    def import_results(
        data_dir: "Path | str",
        hdf5_path: "Path | str",
        file_code: str = "",
        *,
        keep_on_failure: bool = False,
        keep_tmpdir: bool = False,
    ) -> None:
        """Convert Fortran output in *data_dir* and write it to *hdf5_path*.

        Reads the Fortran microscale output (version
        :data:`MICRO_FORTRAN_DATASPEC_VERSION`) from *data_dir*, converts it
        to the HDF5 v2.0.0 format, and writes the resulting datasets into the
        existing HDF5 file at *hdf5_path*.

        The target HDF5 file should already contain ``micro_params`` and
        empty microscale_out datasets (i.e. :meth:`DataStore.create` should
        have been called beforehand).  If a dataset is missing it will be
        created; if it already exists (zero-length from ``DataStore.create``)
        it will be resized and filled.

        Cleanup behaviour:

        - On **success**: *data_dir* is removed unless ``keep_tmpdir=True``.
        - On **failure**: *data_dir* is preserved when
          ``keep_on_failure=True`` *or* ``keep_tmpdir=True``; otherwise it is
          removed.

        :param data_dir: Directory containing the Fortran binary output files
            (``lysis{code}.dat``, ``micro{code}.txt``, ``params.json``, …).
            Typically ``{work_dir}/data/{run_code}/``.
        :type data_dir: Path or str
        :param hdf5_path: Full path to the target ``.h5`` file.
        :type hdf5_path: Path or str
        :param file_code: Output file code suffix used when running the binary
            (e.g. ``"_PLG2_tPA01_Q4"``), defaults to ``""``.
        :type file_code: str, optional
        :param keep_on_failure: Preserve *data_dir* if import raises an
            exception, defaults to ``False``.
        :type keep_on_failure: bool, optional
        :param keep_tmpdir: Always preserve *data_dir*, defaults to ``False``.
        :type keep_tmpdir: bool, optional
        :raises Exception: Re-raises any exception from the read/convert/write
            pipeline after applying the cleanup policy.
        """
        from ..dataio.dataspec import dataspec
        from ..dataio.fileops import read_data_collection
        from ..dataio.dataconvert import convert_data

        data_dir = Path(data_dir)
        hdf5_path = Path(hdf5_path)

        try:
            src_spec = dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"]
            raw = read_data_collection(str(data_dir), [src_spec], [file_code])
            converted = convert_data(raw, MICRO_FORTRAN_DATASPEC_VERSION, "v2.0.0")

            dst_spec = dataspec["v2.0.0"]["microscale_out"]
            with h5py.File(str(hdf5_path), "a") as f:
                for name, ds_spec in dst_spec.data.items():
                    if name not in converted:
                        continue
                    arr = np.asarray(converted[name], dtype=ds_spec.dtype)
                    loc = ds_spec.data_location
                    if loc in f:
                        # Resize existing empty dataset created by DataStore.create()
                        f[loc].resize(arr.shape)
                        f[loc][...] = arr
                    else:
                        # Dataset absent — create it
                        maxshape = tuple(
                            None if (isinstance(s, str) or s < 0) else s
                            for s in ds_spec.shape
                        )
                        f.create_dataset(
                            loc, data=arr, maxshape=maxshape,
                            compression="gzip", dtype=ds_spec.dtype,
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
        hdf5_path: "Path | str",
        *,
        keep_tmpdir: bool = False,
    ) -> None:
        """Execute the complete HDF5-integrated workflow.

        Creates a uniquely named temporary working directory (via
        :func:`tempfile.mkdtemp` with an explicit *dir* argument so that
        ``/tmp`` — which is a ramdisk on many HPC nodes — is never used),
        runs the Fortran binary, imports the results into *hdf5_path*, and
        cleans up.

        The temporary directory is placed alongside the HDF5 file
        (``hdf5_path.parent``) to avoid using ramdisk storage.

        Cleanup policy:

        - **Success**: temporary directory removed unless ``keep_tmpdir=True``.
        - **Failure**: temporary directory is **preserved** for debugging
          (regardless of ``keep_tmpdir``) unless ``keep_tmpdir=False`` AND
          the caller explicitly wants no preservation.  Concretely:
          ``import_results`` is called with ``keep_on_failure=True`` so the
          data survives an import error; the outer ``except`` block removes
          it only when ``keep_tmpdir=False``.

        :param hdf5_path: Full path to the target ``.h5`` file.
        :type hdf5_path: Path or str
        :param keep_tmpdir: Always preserve the temporary directory,
            defaults to ``False``.
        :type keep_tmpdir: bool, optional
        :raises subprocess.CalledProcessError: If the Fortran binary fails.
        :raises Exception: Re-raises any exception from the import pipeline.
        """
        hdf5_path = Path(hdf5_path)
        tmpdir = Path(tempfile.mkdtemp(
            prefix=f"lysis-micro-{self.run.run_code}-",
            dir=str(hdf5_path.parent),
        ))
        try:
            data_dir = self.exec_in_workdir(tmpdir)
            self.import_results(
                data_dir,
                hdf5_path,
                file_code=self.out_file_code,
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

    # def compile(self):
    #     """Compile the Fortran microscale source code (DEPRECATED).
    #
    #     This method attempted to compile the Fortran source using the Intel
    #     Fortran compiler (ifort). It has been disabled because:
    #
    #     - Compilation should be handled separately from execution
    #     - Module loading ("module load oneapi/compiler") doesn't work in subprocess
    #     - Build systems should use proper build tools (Make, CMake, etc.)
    #
    #     Kept as reference for the compilation flags and dependencies.
    #
    #     Original behavior:
    #         - Loaded Intel OneAPI compiler module
    #         - Compiled with: ifort -r8 -mcmodel medium -traceback
    #         - Linked with kiss.o (KISS random number generator)
    #     """
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
