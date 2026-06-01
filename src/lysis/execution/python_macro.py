"""In-process Python (NumPy) macroscale simulation runner.

Provides :class:`PythonRunner` — the in-process counterpart to
:class:`~lysis.execution.fortran.FortranRunner` — and :class:`PythonMacro`,
which runs :class:`~lysis.macroscale.MacroscaleSim` against an HDF5 run file
with no Fortran toolchain required.

Execution is currently **serial**: each of the run's ``macro_simulations``
simulations runs in turn through a single open
:class:`~lysis.dataio.datastore.DataStore`, because HDF5 cannot have two
simulations writing the same file concurrently.  Per-simulation seeds are still
derived with :class:`numpy.random.SeedSequence` (matching
:meth:`~lysis.execution.fortran_macro.FortranMacro.exec_in_workdir`), so a
future dispatch-collect / MPI4Py parallel implementation will produce identical
results.  See GitHub issues #59 (dispatch-collect) and #60 (MPI4Py).
"""

import dataclasses

from abc import abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..config.constants import CONST
from ..config.run import Run
from ._state import require_macro_empty
from .base import SimulationRunner


@dataclass
class PythonRunner(SimulationRunner):
    """Abstract base for in-process (pure-Python) simulation runners.

    The Python counterpart to :class:`~lysis.execution.fortran.FortranRunner`:
    instead of launching a compiled binary it drives an in-process simulation
    class (e.g. :class:`~lysis.macroscale.MacroscaleSim`) that writes its
    results straight into the run's HDF5 file.

    Subclasses supply the scale-specific details via the abstract hooks.

    :ivar run: The Run object containing all simulation parameters and data.
    :vartype run: Run
    """

    run: Run = None

    # ------------------------------------------------------------------
    # Abstract hooks (implemented by subclasses)
    # ------------------------------------------------------------------

    @abstractmethod
    def _simulation_class(self):
        """Return the in-process simulation class to instantiate per simulation.

        The class must accept ``(run, sim_number=...)`` and expose a ``go()``
        method that writes results into ``run.data``.

        :rtype: type
        """
        ...

    @abstractmethod
    def _get_params(self):
        """Return the relevant parameters instance from the Run.

        :return: The parameters object (e.g. ``self.run.macro_params``).
        """
        ...

    @abstractmethod
    def _set_params(self, params) -> None:
        """Assign the relevant parameters instance back onto the Run.

        :param params: The parameters object to store on the Run.
        """
        ...

    @abstractmethod
    def _seed_field(self) -> str:
        """Return the name of the seed field on the parameters object.

        :rtype: str
        """
        ...

    @abstractmethod
    def _simulations_field(self) -> str:
        """Return the name of the simulation-count field on the parameters.

        :rtype: str
        """
        ...

    @abstractmethod
    def _scale(self) -> str:
        """Return the provenance scale name (``"micro"`` or ``"macro"``).

        :rtype: str
        """
        ...

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_hdf5(cls, hdf5_path: "Path | str") -> "PythonRunner":
        """Construct a runner from an existing HDF5 run file.

        Validates the file is in the ``MACRO_EMPTY`` state (via
        :func:`~lysis.execution._state.require_macro_empty`), then rebuilds the
        :class:`~lysis.config.run.Run` and loads its parameters from the file.

        :param hdf5_path: Full path to the ``.h5`` file (must contain both
            ``micro_params`` and ``macro_params``).
        :type hdf5_path: Path or str
        :return: A fully configured runner instance.
        :rtype: PythonRunner
        :raises ValueError: If the file is not in the ``MACRO_EMPTY`` state.
        """
        hdf5_path = Path(hdf5_path)
        require_macro_empty(hdf5_path)
        run = Run(str(hdf5_path.parent), run_code=hdf5_path.stem)
        run.load_params_from_hdf5()
        return cls(run=run)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def execute(self) -> None:
        """Run every simulation for this run, in series, into its HDF5 file.

        Opens the run's :class:`~lysis.dataio.datastore.DataStore` in append
        mode (a single writer), then for each of the ``macro_simulations``
        simulations:

        1. Derives an independent seed from
           ``np.random.SeedSequence(<base seed>).generate_state(n)``.
        2. Replaces the (frozen) parameters object with one carrying that seed.
        3. Runs the in-process simulation, which writes its own
           per-simulation slot via the simulation's ``go()`` /
           ``record_data_to_disk()``.

        Finally stamps provenance onto the per-scale params group -- the
        ``pipeline_*`` family plus the ``backend_*`` family with
        ``backend_type="python"`` (#77) -- then closes the file.

        Serial execution is a deliberate interim measure (HDF5 single-writer
        limitation); see issues #59 / #60.
        """
        orig_params = self._get_params()
        seed_field = self._seed_field()
        n_sims = getattr(orig_params, self._simulations_field())
        # The persisted seed is a canonical entropy string (#97/#109); decode
        # it to an int before splitting, matching the Fortran macro path.
        seeds = np.random.SeedSequence(
            orig_params.seed_as_int()
        ).generate_state(n_sims)

        sim_class = self._simulation_class()
        ds = self.run.open_data(mode="a")
        try:
            for sim_number in range(n_sims):
                self._set_params(
                    dataclasses.replace(
                        orig_params, **{seed_field: int(seeds[sim_number])}
                    )
                )
                sim_class(self.run, sim_number=sim_number).go()

            # Restore the original (unmutated) parameters object.
            self._set_params(orig_params)

            # Stamp the v1.0.0 provenance schema: the pipeline_* family (the
            # src/lysis wrapper) plus the backend_* family with
            # backend_type="python" (#77).  Mirrors what import_collection
            # stamps for the Fortran backend.
            from ..tools.provenance import gather_python_backend_provenance
            ds.stamp_provenance(self._scale(), "pipeline")
            ds.stamp_provenance(
                self._scale(),
                "backend",
                replace_backend_attrs=gather_python_backend_provenance(),
                backend_type="python",
            )
            ds._file.flush()
        finally:
            ds.close()
            self.run.data = None


class PythonMacro(PythonRunner):
    """In-process NumPy macroscale runner.

    Drives :class:`~lysis.macroscale.MacroscaleSim` against an HDF5 run file.
    """

    def _simulation_class(self):
        from ..macroscale import MacroscaleSim  # noqa: PLC0415
        return MacroscaleSim

    def _get_params(self):
        return self.run.macro_params

    def _set_params(self, params) -> None:
        self.run.macro_params = params

    def _seed_field(self) -> str:
        return "macro_seed"

    def _simulations_field(self) -> str:
        return "macro_simulations"

    def _scale(self) -> str:
        return "macro"
