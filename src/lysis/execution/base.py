"""Abstract base class for all simulation runners.

Defines :class:`SimulationRunner`, the common interface shared by every
class that executes a fibrinolysis simulation — whether via a compiled
Fortran subprocess, an in-process NumPy implementation, or future GPU-
accelerated code.
"""

from abc import ABC, abstractmethod

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class SimulationRunner(ABC):
    """Abstract base class for all simulation execution classes.

    Defines the common interface that all simulation runners must implement,
    regardless of whether they execute Fortran subprocesses, run in-process
    Python simulations, or dispatch to GPU-accelerated code.

    Subclasses must implement :meth:`execute` to perform the actual
    simulation work.
    """

    @abstractmethod
    def execute(self) -> None:
        """Execute the simulation.

        Implementations should run the simulation to completion, writing
        results to the appropriate output location (HDF5 file, working
        directory, etc.).
        """
        ...
