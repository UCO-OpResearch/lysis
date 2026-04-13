"""Wrapper for executing the Fortran macroscale simulation binary.

Provides :class:`FortranMacro`, which manages execution of the compiled
Fortran macroscale simulation.  It converts Python
:class:`~lysis.config.run.Run` parameters to command-line arguments,
generates the required neighbor-structure input file, and captures binary
stdout to a log file.
"""

import os

from dataclasses import dataclass
from typing import AnyStr

from ..config.parameters import MacroParameters, MicroParameters
from ..geometry.edge_grid import generate_fortran_neighborhood_structure
from .fortran import FortranRunner

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
    command-line arguments, generates the required neighbor-structure input
    file (``neighbors.dat``), and captures binary stdout to a log file.

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

    # ------------------------------------------------------------------
    # Neighborhood generation
    # ------------------------------------------------------------------

    def generate_neighborhoods(self):
        """Generate and write the neighbor-structure file for Fortran.

        Pre-calculates the 8-neighbor structure for every edge in the
        hexagonal grid and writes it to ``neighbors.dat`` in
        ``run.os_path``.  Fortran uses 1-based indexing, so all indices are
        incremented by 1.

        :raises OSError: If the output file cannot be written.
        """
        fort_neighbors = (
            generate_fortran_neighborhood_structure(
                self.run.macro_params.rows, self.run.macro_params.cols
            )
            + 1  # Convert 0-based (Python) → 1-based (Fortran)
        )
        fort_neighbors.tofile(
            os.path.join(self.run.os_path, "neighbors.dat"), sep=os.linesep
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def execute(self) -> None:
        """Execute the Fortran macroscale simulation.

        Performs the complete workflow:

        1. Generates ``neighbors.dat`` via :meth:`generate_neighborhoods`.
        2. Builds command-line arguments via :meth:`exec_command`.
        3. Runs the Fortran binary, capturing stdout to
           ``macro{out_file_code}.txt`` in ``run.os_path``.

        Blocks until the Fortran binary exits.

        :raises FileNotFoundError: If the executable does not exist.
        :raises OSError: If output files cannot be created.

        Side effects:
            - Writes ``neighbors.dat`` to ``run.os_path``.
            - Writes ``macro{out_file_code}.txt`` to ``run.os_path``.
        """
        self.generate_neighborhoods()
        super().execute()
