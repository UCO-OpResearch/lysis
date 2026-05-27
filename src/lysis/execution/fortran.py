"""Abstract base class for Fortran subprocess simulation runners.

Provides :class:`FortranRunner`, which captures the logic shared by all
classes that execute a compiled Fortran simulation binary:

- Converting Python :class:`~lysis.config.parameters.Parameters` to
  command-line arguments (handling units, index offsets, scaling, and
  RNG seed splitting for parallel runs).
- Launching the Fortran binary via :mod:`subprocess` and capturing its
  stdout to a log file.
- Importing Fortran output into the run's HDF5 file.
- Running the complete HDF5-integrated workflow.

Concrete subclasses (:class:`~lysis.execution.fortran_macro.FortranMacro`
and :class:`~lysis.execution.fortran_micro.FortranMicro`) provide the
parameter-class-specific hooks.
"""

import inspect
import os
import shutil
import subprocess
import tempfile
import warnings

from abc import abstractmethod
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import AnyStr

import numpy as np

from pint import Quantity

from ..config.run import Run
from ..dataio.datastore import DataStore
from .base import SimulationRunner


#: Dataspec version produced by the Fortran microscale binary.
#: Update this constant when the binary's output format changes.
MICRO_FORTRAN_DATASPEC_VERSION: str = "v1.99.0"

#: Dataspec version produced by the Fortran macroscale binary.
#: Update this constant when the binary's output format changes.
MACRO_FORTRAN_DATASPEC_VERSION: str = "v1.99.0"


@dataclass
class FortranRunner(SimulationRunner):
    """Abstract base for classes that execute compiled Fortran simulation binaries.

    Encapsulates the logic shared by :class:`~lysis.execution.fortran_macro.FortranMacro`
    and :class:`~lysis.execution.fortran_micro.FortranMicro`:

    - Converting :class:`~lysis.config.parameters.Parameters` fields to
      command-line arguments via :meth:`exec_command`.
    - Running the subprocess and capturing stdout via :meth:`execute`.

    Subclasses must implement all abstract hook methods to supply the
    parameter-class-specific details used by the template methods above.

    :ivar run: The Run object containing all simulation parameters and data.
    :vartype run: Run
    :ivar cwd: Working directory for executing the Fortran binary.
    :vartype cwd: str
    :ivar executable: Path to the compiled Fortran binary.
    :vartype executable: str
    :ivar out_file_code: Code suffix for output files (e.g. ``_PLG2_tPA01_Q4``).
    :vartype out_file_code: str
    :ivar index: Optional index for parallel runs.  When set, the RNG seed is
        split so each index produces an independent stream and
        ``__{index:02}`` is appended to ``out_file_code``.  The simulation
        count for this task is determined by :attr:`num_children`: when
        ``num_children`` is ``None`` (legacy single-sim-per-task path) the
        count is forced to 1; otherwise it is partitioned across
        ``num_children`` tasks (see :meth:`exec_command`).
    :vartype index: int or None
    :ivar num_children: Optional number of sibling tasks sharing this run.
        When set with :attr:`index`, the parameter ``{simulations}`` is
        partitioned across ``num_children`` tasks (with the remainder
        distributed across the lowest-indexed tasks).  Subclasses must
        return ``num_children`` from :meth:`_seed_split_count` so all
        sibling tasks draw seeds from a stream of identical size.  Defaults
        to ``None`` (preserves legacy single-sim-per-task behavior).
    :vartype num_children: int or None
    """

    run: Run = None
    cwd: AnyStr = "."
    executable: AnyStr = None
    out_file_code: AnyStr = ""
    index: int = None
    num_children: int = None
    #: When True, a mismatch between the Fortran binary's embedded build
    #: stamp and the current ``src/fortran/`` source state downgrades from
    #: a :class:`~lysis.tools.provenance.StaleBinaryError` to a loud
    #: :class:`UserWarning` plus a banner written to each Fortran stdout
    #: log file and override metadata stamped onto the resulting HDF5
    #: group.  Defaults to False; :meth:`_verify_binary_version` also
    #: honours the ``LYSIS_ALLOW_STALE_BINARY`` env var so callers don't
    #: need to combine the two override sources themselves.
    allow_stale_binary: bool = False
    #: When True, :meth:`_verify_binary_version` is bypassed entirely —
    #: no ``<binary> --version`` subprocess, no staleness check, no
    #: warning banner.  Set by the historical-build workflow
    #: (``lysis run-{micro,macro} --fortran-commit ...``) where the
    #: binary↔``src/fortran/`` mismatch is intentional and the historical
    #: binary may not implement ``--version`` at all.  Pair with
    #: :attr:`historical_binary_attrs` so the synthesised binary
    #: provenance is still stamped into the HDF5 file.
    skip_binary_verification: bool = False
    #: Optional pre-computed binary-provenance dict shipped to
    #: :meth:`~lysis.dataio.datastore.DataStore.import_collection` as
    #: ``replace_binary_attrs``.  Populated by the historical-build
    #: workflow from
    #: :func:`~lysis.tools.provenance.gather_historical_binary_provenance`.
    #: When ``None``, the default
    #: :func:`~lysis.tools.provenance.gather_binary_provenance` path is
    #: used (queries ``<binary> --version``).
    historical_binary_attrs: "dict | None" = None

    # ------------------------------------------------------------------
    # Abstract hooks (implemented by subclasses)
    # ------------------------------------------------------------------

    @abstractmethod
    def _get_params(self):
        """Return the relevant parameters instance from the Run.

        :return: The parameters object (e.g. ``self.run.micro_params``).
        """
        ...

    @abstractmethod
    def _get_params_class(self):
        """Return the Parameters class used by this runner.

        :return: The parameters class (e.g. :class:`~lysis.config.parameters.MicroParameters`).
        :rtype: type
        """
        ...

    @abstractmethod
    def _seed_field(self) -> str:
        """Return the name of the seed field in the parameters dict.

        :return: Field name, e.g. ``"micro_seed"`` or ``"macro_seed"``.
        :rtype: str
        """
        ...

    @abstractmethod
    def _simulations_field(self) -> str:
        """Return the name of the simulation-count field in the parameters dict.

        :return: Field name, e.g. ``"micro_simulations"`` or ``"macro_simulations"``.
        :rtype: str
        """
        ...

    @abstractmethod
    def _seed_split_count(self, params: dict) -> int:
        """Return the count for :meth:`numpy.random.SeedSequence.generate_state`.

        Called *before* ``params`` is mutated (i.e. before the simulations
        field is partitioned), so implementations may safely read ``params``.

        Implementations **must** return a value that is identical across
        sibling tasks sharing the same base seed; otherwise sibling tasks
        will draw seeds from streams of different sizes and the partition
        will not be reproducibility-correct.  The standard implementation
        returns :attr:`num_children` when set (array path) and falls back
        to ``self.index + 1`` for the legacy single-sim-per-task path,
        which preserves bit-for-bit reproducibility of pre-existing
        single-task artifacts.

        :param params: The parameters dict (from :func:`dataclasses.asdict`).
        :type params: dict
        :return: The number of seeds to generate.
        :rtype: int
        """
        ...

    @abstractmethod
    def _base_arguments(self) -> list:
        """Return the opening CLI arguments before the parameter loop.

        Typically ``["--runCode", run_code, "--outFileCode", out_file_code]``
        plus any additional fixed arguments (e.g. ``--inFileCode`` for macro).

        :return: List of CLI argument strings.
        :rtype: list[str]
        """
        ...

    def _post_arguments(self, params: dict) -> list:
        """Return any trailing CLI arguments appended after the parameter loop.

        The default implementation returns an empty list.  Subclasses may
        override to inject extra arguments (e.g. ``FortranMacro`` appends
        ``--radius`` and ``--bs`` from ``micro_params``).

        :param params: The parameters dict (from :func:`dataclasses.asdict`).
        :type params: dict
        :return: List of CLI argument strings (may be empty).
        :rtype: list[str]
        """
        return []

    @abstractmethod
    def _log_prefix(self) -> str:
        """Return the prefix for the stdout log filename.

        The log file is written to ``{run.os_path}/{prefix}{out_file_code}.txt``.

        :return: Prefix string, e.g. ``"micro"`` or ``"macro"``.
        :rtype: str
        """
        ...

    @abstractmethod
    def _collection_name(self) -> str:
        """Return the DataStore collection name for this runner's output.

        :return: Collection name, e.g. ``"microscale_out"`` or
            ``"macroscale_out"``.
        :rtype: str
        """
        ...

    @abstractmethod
    def _fortran_dataspec_version(self) -> str:
        """Return the dataspec version produced by this runner's Fortran binary.

        :return: Version string, e.g. :data:`MICRO_FORTRAN_DATASPEC_VERSION`.
        :rtype: str
        """
        ...

    def _pre_execute(self) -> None:
        """Hook called by :meth:`execute` before running the subprocess.

        The default implementation is a no-op.  Subclasses may override
        to perform setup required by the legacy ``execute()`` path (e.g.
        generating input files).
        """

    def _verify_binary_version(self) -> dict:
        """Preflight-check that :attr:`executable` matches ``src/fortran/``.

        Delegates to
        :func:`~lysis.tools.provenance.verify_binary_matches_source`
        and caches the result so multi-simulation paths (e.g. the macro
        loop) issue the warning exactly once even when called per-sim.

        When :attr:`skip_binary_verification` is True (historical-build
        workflow), this is a no-op — no subprocess is spawned, no
        warning banner is generated.  The caller has independently
        synthesised the binary provenance and will pass it through
        :attr:`historical_binary_attrs`.

        :return: The verify dict — empty on match or when skipped; on
            overridden mismatch contains ``"banner"`` (text to prepend
            to the Fortran stdout log) plus the ``stale_binary_override``
            flag forwarded to
            :meth:`~lysis.dataio.datastore.DataStore.import_collection`
            via :meth:`_binary_hdf5_attrs`.  The binary's
            commit/dirty/compiler stamps are no longer included here —
            they are gathered fresh and stamped unconditionally by
            :meth:`~lysis.dataio.datastore.DataStore.stamp_provenance`
            on every run.
        :rtype: dict
        :raises ~lysis.tools.provenance.StaleBinaryError: When the
            binary stamp disagrees with ``src/fortran/`` and
            :attr:`allow_stale_binary` is False.
        """
        if self.skip_binary_verification:
            self._binary_check_info = {}
            return self._binary_check_info
        if getattr(self, "_binary_check_info", None) is None:
            from ..tools.provenance import (  # noqa: PLC0415
                allow_stale_from_env,
                verify_binary_matches_source,
            )

            allow_stale = self.allow_stale_binary or allow_stale_from_env()
            self._binary_check_info = verify_binary_matches_source(
                self.executable,
                allow_stale=allow_stale,
            )
        return self._binary_check_info

    @property
    def _binary_hdf5_attrs(self) -> dict:
        """Override-only attrs from :meth:`_verify_binary_version` to merge
        into the binary stamp; ``{}`` when the binary matched the source.

        Today this carries only :data:`CONST.STALE_BINARY_OVERRIDE_ATTR`
        when the staleness check was overridden.  The binary's commit /
        dirty / compiler stamps are written unconditionally by
        :meth:`~lysis.dataio.datastore.DataStore.stamp_provenance` —
        they are no longer routed through this property.
        """
        info = getattr(self, "_binary_check_info", None) or {}
        return {k: v for k, v in info.items() if k != "banner"}

    def _write_stale_banner(self, fh) -> None:
        """Write the stale-binary banner to *fh* if the preflight was overridden.

        No-op on a clean match.  Must be called on the same file handle
        the subprocess will write its stdout to, BEFORE
        :func:`subprocess.run`, so the banner lands at the top of the log.

        :param fh: An open writable text file handle.
        """
        info = getattr(self, "_binary_check_info", None) or {}
        banner = info.get("banner")
        if banner:
            fh.write(banner)
            fh.flush()

    # ------------------------------------------------------------------
    # Template methods
    # ------------------------------------------------------------------

    def _params_to_arguments(self, params: dict) -> list:
        """Convert a parameters dict to a flat list of CLI argument pairs.

        Handles unit conversions for :class:`pint.Quantity` values, index
        adjustments (``"-1"`` Fortran-name suffix: Python 0-based →
        Fortran 1-based), scaling (``"*100"`` suffix), and skips parameters
        that equal their class default.

        :param params: Parameter dict produced by :func:`dataclasses.asdict`
            from the relevant parameters dataclass.  The seed field should
            be a ``np.uint32``; it is reinterpreted to a signed ``INTEGER*4``
            literal here via the ``|uint32`` Fortran-name suffix.
        :type params: dict
        :return: Flat list of CLI argument pairs, e.g.
            ``["--rows", "10", "--cols", "20"]``.
        :rtype: list[str]
        """
        arguments = []
        params_class = self._get_params_class()
        sig = inspect.signature(params_class)
        fortran_names = params_class.fortran_names()
        units = params_class.units()

        for key in sig.parameters:
            if key in fortran_names and params[key] != sig.parameters[key].default:
                value = params[key]
                if isinstance(value, Quantity):
                    value = value.m_as(units[key])

                fname = fortran_names[key]
                if fname.endswith("|uint32"):
                    # Fortran read(*, *) expects a signed INTEGER*4 literal.
                    # Reinterpret the Python np.uint32 bits as np.int32
                    # (.astype is the documented, unambiguous bit-preserving
                    # cast for same-width dtype changes; constructor coercion
                    # can raise OverflowError in NumPy 2.x). The C KISS side
                    # then reads the same bytes back as uint_least32_t, so a
                    # high-bit-set uint32 round-trips exactly.
                    signed = int(
                        np.array(value, dtype=np.uint32).astype(np.int32)
                    )
                    arguments += ["--" + fname[:-7], str(signed)]
                elif fname.endswith("-1"):
                    arguments += ["--" + fname[:-2], str(value + 1)]
                elif fname.endswith("*100"):
                    arguments += ["--" + fname[:-4], str(value // 100)]
                else:
                    arguments += ["--" + fname, str(value)]

        return arguments

    def exec_command(self) -> list:
        """Build the command-line arguments for the Fortran binary.

        Converts the relevant parameter object to a list of command-line
        arguments compatible with the Fortran program.  Handles:

        - RNG seed splitting when :attr:`index` is set.
        - Unit conversions for :class:`pint.Quantity` objects.
        - Index adjustments for parameters with a ``"-1"`` Fortran-name suffix
          (Python 0-based → Fortran 1-based).
        - Scaling for parameters with a ``"*100"`` Fortran-name suffix.
        - Only non-default parameters are included.

        .. note::
            When :attr:`index` is not ``None`` this method mutates
            :attr:`out_file_code` by appending ``__{index:02}``.  Call it
            exactly once per execution.

        :return: Complete command as ``[executable, arg1, val1, ...]``
            suitable for :func:`subprocess.run`.
        :rtype: list[str]
        """
        params = asdict(self._get_params())

        # Seed splitting for parallel (multi-node) runs
        if self.index is not None:
            stream = np.random.SeedSequence(params[self._seed_field()])
            split_count = self._seed_split_count(params)  # read before mutation
            seeds = stream.generate_state(split_count)
            sim_field = self._simulations_field()
            if self.num_children is None:
                # Legacy: one simulation per task.
                params[sim_field] = 1
            else:
                # Array partition: distribute total sims across children,
                # giving the lowest-indexed tasks the remainder so
                # sum(per_task) == total.
                total = int(params[sim_field])
                k = self.num_children
                params[sim_field] = total // k + (1 if self.index < total % k else 0)
            params[self._seed_field()] = seeds[self.index]
            self.out_file_code = self.out_file_code + f"__{self.index:02}"

        arguments = self._base_arguments()
        arguments += self._params_to_arguments(params)
        arguments += self._post_arguments(params)
        return [self.executable] + arguments

    def execute(self) -> None:
        """Execute the Fortran binary and capture its stdout to a log file.

        .. deprecated::
            Use :meth:`exec_in_workdir` / :meth:`run_full` for the
            HDF5-integrated workflow instead.

        Builds the command via :meth:`exec_command`, then runs the binary
        with ``cwd=self.cwd``, redirecting stdout to
        ``{run.os_path}/{prefix}{out_file_code}.txt``.

        Blocks until the Fortran binary exits.

        :raises FileNotFoundError: If the executable does not exist.
        :raises subprocess.CalledProcessError: If the binary exits with a
            non-zero status (only when :func:`subprocess.run` is called with
            ``check=True``; the current implementation does not set it).
        :raises OSError: If the log file cannot be created.
        """
        warnings.warn(
            f"{type(self).__name__}.execute() is deprecated. "
            "Use exec_in_workdir() / run_full() for the HDF5-integrated "
            "workflow instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._pre_execute()
        command = self.exec_command()
        output_file_name = os.path.join(
            self.run.os_path,
            self._log_prefix() + self.out_file_code + ".txt",
        )
        with open(output_file_name, "w") as file:
            subprocess.run(command, stdout=file, cwd=self.cwd)

    # ------------------------------------------------------------------
    # Import and full workflow
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
        the Fortran output (version :meth:`_fortran_dataspec_version`),
        convert it to HDF5 v2.0.0, and populate the empty collection datasets
        (named by :meth:`_collection_name`).

        The target HDF5 file should already contain the collection with empty
        datasets (i.e. :meth:`~lysis.dataio.datastore.DataStore.create` or
        :meth:`~lysis.dataio.datastore.DataStore.initialize_macroscale` should
        have been called beforehand).

        Cleanup behaviour:

        - On **success**: *data_dir* is removed unless ``keep_tmpdir=True``.
        - On **failure**: *data_dir* is preserved when
          ``keep_on_failure=True`` *or* ``keep_tmpdir=True``; otherwise it
          is removed.

        :param data_dir: Directory containing the Fortran binary output files.
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
                    self._collection_name(),
                    self._fortran_dataspec_version(),
                    str(data_dir),
                    [self.out_file_code],
                    binary_executable=self.executable,
                    binary_override=self._binary_hdf5_attrs,
                    replace_binary_attrs=self.historical_binary_attrs,
                )

            if not keep_tmpdir:
                shutil.rmtree(data_dir)

        except Exception:
            if not keep_on_failure and not keep_tmpdir:
                shutil.rmtree(data_dir, ignore_errors=True)
            raise

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
                prefix=f"lysis-{self._log_prefix()}-{self.run.run_code}-",
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
