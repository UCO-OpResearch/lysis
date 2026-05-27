"""Shared HDF5 lifecycle-state precondition for macroscale runners.

Both :class:`~lysis.execution.fortran_macro.FortranMacro` and
:class:`~lysis.execution.python_macro.PythonMacro` require the run's HDF5 file
to be in the :attr:`~lysis.dataio.datastore.HDF5State.MACRO_EMPTY` state before
they run.  Centralising the check here keeps the two backends' error messages
identical.
"""

from pathlib import Path

from ..dataio.datastore import DataStore, HDF5State


def require_macro_empty(hdf5_path: "Path | str") -> None:
    """Confirm *hdf5_path* is ready for a macroscale run.

    Opens the file read-only and verifies its lifecycle state is
    :attr:`~lysis.dataio.datastore.HDF5State.MACRO_EMPTY` — microscale output
    present and macroscale initialised but not yet run.

    :param hdf5_path: Path to the run's ``.h5`` file.
    :type hdf5_path: Path or str
    :raises ValueError: If microscale data is missing or inconsistent, if
        macroscale has not been initialised, or if macroscale results already
        exist.
    """
    hdf5_path = Path(hdf5_path)
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
