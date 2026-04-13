"""Abstract base class for Fortran subprocess simulation runners.

Provides :class:`FortranRunner`, which captures the logic shared by all
classes that execute a compiled Fortran simulation binary:

- Converting Python :class:`~lysis.config.parameters.Parameters` to
  command-line arguments (handling units, index offsets, scaling, and
  RNG seed splitting for parallel runs).
- Launching the Fortran binary via :mod:`subprocess` and capturing its
  stdout to a log file.

Concrete subclasses (:class:`~lysis.execution.fortran_macro.FortranMacro`
and :class:`~lysis.execution.fortran_micro.FortranMicro`) provide the
parameter-class-specific hooks.
"""

import inspect
import os
import subprocess

from abc import abstractmethod
from dataclasses import asdict, dataclass
from typing import AnyStr

import numpy as np

from pint import Quantity

from ..config.run import Run
from .base import SimulationRunner

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


#: Dataspec version produced by the Fortran microscale binary.
#: Update this constant when the binary's output format changes.
MICRO_FORTRAN_DATASPEC_VERSION: str = "v1.99.0"


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
        split so each index produces an independent stream, ``{simulations}``
        is forced to 1, and ``__{index:02}`` is appended to ``out_file_code``.
    :vartype index: int or None
    """

    run: Run = None
    cwd: AnyStr = "."
    executable: AnyStr = None
    out_file_code: AnyStr = ""
    index: int = None

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
        field is set to 1), so implementations may safely read ``params``.

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

    # ------------------------------------------------------------------
    # Template methods
    # ------------------------------------------------------------------

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
            params[self._simulations_field()] = 1
            params[self._seed_field()] = seeds[self.index]
            self.out_file_code = self.out_file_code + f"__{self.index:02}"

        # Cast seed to int32 to match Fortran INTEGER*4
        params[self._seed_field()] = int(
            np.array(params[self._seed_field()]).astype(np.int32)
        )

        arguments = self._base_arguments()

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
                if fname.endswith("-1"):
                    arguments += ["--" + fname[:-2], str(value + 1)]
                elif fname.endswith("*100"):
                    arguments += ["--" + fname[:-4], str(value // 100)]
                else:
                    arguments += ["--" + fname, str(value)]

        arguments += self._post_arguments(params)
        return [self.executable] + arguments

    def execute(self) -> None:
        """Execute the Fortran binary and capture its stdout to a log file.

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
        command = self.exec_command()
        output_file_name = os.path.join(
            self.run.os_path,
            self._log_prefix() + self.out_file_code + ".txt",
        )
        with open(output_file_name, "w") as file:
            subprocess.run(command, stdout=file, cwd=self.cwd)
