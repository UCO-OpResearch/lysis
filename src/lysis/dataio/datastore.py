"""Read-only HDF5 interface for fibrinolysis simulation data (v2.0.0 dataspec).

This module provides the :class:`DataStore` class, the primary entry point for
reading and writing simulation output files that conform to the lysis HDF5 v2.0.0
data specification.  It also exposes supporting classes that give structured,
dot-access navigation of the collections and datasets stored inside those files.

Data model
----------
An HDF5 file managed by this module may contain up to three data collections,
each corresponding to a stage of the simulation pipeline:

- **microscale_out** — Output from the microscale (C) model.  All simulations
  are stored in a single combined dataset group (``simulations_combined=True``).
  Access datasets directly via dot-notation::

      ds.microscale_out.tpa_leaving_time[:]

- **macroscale_out** — Output from the macroscale (Fortran) model.  One dataset
  group per simulation (``simulations_combined=False``).  Index by simulation
  number to get a :class:`SimulationView`, then access datasets on it::

      ds.macroscale_out[3].fiber_degrade_time[:]

- **macroscale_in** — Derived input for the macroscale model, computed lazily
  from ``microscale_out`` data via
  :func:`~lysis.dataio.dataconvert.generate_macroscale_in`.
  Stored in memory (not on disk) and generated on first access.

Parameters (``MicroParameters``, ``MacroParameters``) are stored as HDF5
attributes and loaded automatically when the corresponding collection is present.

Classes
-------
DataStore
    Primary interface.  Opens, validates, and navigates an HDF5 file.
    Supports the context-manager protocol (``with DataStore(...) as ds:``).

DataCollection
    Lazy-access wrapper around one on-disk collection.  Returns
    ``h5py.Dataset`` objects (combined) or :class:`SimulationView` objects
    (per-simulation).

SimulationView
    Dot-access to the datasets belonging to a single simulation within a
    per-simulation collection.

DerivedDataCollection
    In-memory wrapper for computed collections.  Identical dot-access
    interface to :class:`DataCollection` but backed by numpy arrays.

DataStatus
    :class:`~enum.Flag` tracking the lifecycle state of a DataStore
    (``INITIALIZED``, ``LOADED``, ``SAVED``, ``FILLED``).

Typical usage
-------------
Open an existing file for reading::

    from lysis.dataio.datastore import DataStore

    with DataStore("run01", "/path/to/data") as ds:
        leaving_times = ds.microscale_out.tpa_leaving_time[:]
        degrade_time  = ds.macroscale_out[0].fiber_degrade_time[:]
        rows          = ds.macro_params.rows

Create a new file and populate it::

    ds = DataStore.create("run01", "/path/to/data", micro_params)
    # … write microscale results …
    ds.initialize_macroscale(macro_params)
    # … write macroscale results …
    ds.close()

Module-level constant
---------------------
COMPATIBLE_DATASPEC_VERSION
    The dataspec version string this module is built against (``"v2.0.0"``).
    A :class:`UserWarning` is raised at import time if the active dataspec
    ``"hdf5"`` tag points to a different version.
"""
import dataclasses
import os
import warnings

from enum import Enum, Flag, auto, unique
from typing import AnyStr

import numpy as np
import h5py

from ..config.constants import CONST
from ..config.parameters import MicroParameters, MacroParameters
from .dataspec import DataCollectionSpec, DataSetSpec, dataspec, parse_shape, tags
from .fileops import read_dataset, write_dataset, _validate_hdf5_version, read_data_collection


#: The dataspec version that this module is compatible with.
COMPATIBLE_DATASPEC_VERSION = "v2.0.0"

if dataspec["hdf5"].version != COMPATIBLE_DATASPEC_VERSION:
    warnings.warn(
        f"DataStore is built for dataspec '{COMPATIBLE_DATASPEC_VERSION}', "
        f"but the 'hdf5' tag currently points to '{dataspec['hdf5'].version}'. "
        f"DataStore may not work correctly with the current HDF5 spec.",
        UserWarning,
        stacklevel=1,
    )


class ImportCollectionError(RuntimeError):
    """Raised when :meth:`DataStore.import_collection` fails during the write
    phase.

    Signals that an import failed *after* HDF5 mutation began and that the
    target collection has been rolled back to its empty, freshly-initialized
    state (``MICRO_EMPTY`` / ``MACRO_EMPTY``).  The original failure is chained
    as the exception's ``__cause__``.
    """


@unique
class DataStatus(Flag):
    NONE = 0
    INITIALIZED = auto()
    LOADED = auto()
    SAVED = auto()
    FILLED = auto()


@unique
class HDF5State(Enum):
    """Lifecycle state of a lysis HDF5 run file.

    Describes which simulation stages have been completed, as inferred from
    which data groups are present and whether their datasets are populated.

    States follow the pipeline order:

    ``MICRO_EMPTY`` → ``MICRO_FILLED`` → ``MACRO_EMPTY`` → ``MACRO_FILLED``

    :cvar MICRO_EMPTY: Microscale group present with empty datasets.
        File is in the post-:meth:`~DataStore.create` (init-experiment) state,
        ready for the microscale simulation to be run.
    :cvar MICRO_FILLED: Microscale datasets populated; no macroscale group.
        File is in the post-run-micro state, ready for
        :meth:`~DataStore.initialize_macroscale`.
    :cvar MACRO_EMPTY: Macroscale group present with empty datasets; microscale
        datasets populated.  File is in the post-:meth:`~DataStore.initialize_macroscale`
        state, ready for the macroscale simulation to be run.
    :cvar MACRO_FILLED: Both microscale and macroscale datasets populated.
        File is in the post-run-macro (complete) state.
    :cvar INCONSISTENT: Macroscale group present but microscale datasets are
        empty.  This state should not arise in normal usage; it indicates the
        HDF5 file has been modified outside the normal pipeline.
    """

    MICRO_EMPTY = auto()
    MICRO_FILLED = auto()
    MACRO_EMPTY = auto()
    MACRO_FILLED = auto()
    INCONSISTENT = auto()


def h5_tree(val: h5py.Dataset, pre: AnyStr = "") -> str:
    """Recursively prints the tree of an HDF5 file's contents.

    Copied from https://stackoverflow.com/questions/61133916/is-there-in-python-a-single-function-that-shows-the-full-structure-of-a-hdf5-fi

    Args:
        val: The item in the HDF5 to print
        pre: The current indentation
    """
    output = ""
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                output += pre + "└── " + key + os.linesep
                output += h5_tree(val, pre + "    ")
            else:
                try:
                    output += pre + "└── " + key + f" {val.shape}" + os.linesep
                except TypeError:
                    output += pre + "└── " + key + " (scalar)" + os.linesep
        else:
            if type(val) == h5py._hl.group.Group:
                output += pre + "├── " + key + os.linesep
                output += h5_tree(val, pre + "│   ")
            else:
                try:
                    output += pre + "├── " + key + f" {val.shape}" + os.linesep
                except TypeError:
                    output += pre + "├── " + key + " (scalar)" + os.linesep
    return output


def _group_path_from_spec(collection_spec):
    """Extract the top-level HDF5 group path from a collection spec.

    Looks at the collection's params ``data_location`` (e.g., ``"micro_data"``)
    to determine the HDF5 group that holds this collection's data.

    :param collection_spec: The collection specification to inspect.
    :type collection_spec: DataCollectionSpec
    :return: The HDF5 group path, or ``None`` if not determinable.
    :rtype: str | None
    """
    if collection_spec.params is not None and collection_spec.params.data_location:
        return collection_spec.params.data_location
    return None


class SimulationView:
    """Structured view of one simulation's datasets within a per-sim collection.

    Provides dot-access to HDF5 datasets for a single simulation index.
    Dataset names are resolved via the spec's ``data_location`` field,
    formatted with the simulation index.  Whether datasets are writable
    depends on the mode of the underlying :class:`h5py.File`; this class
    imposes no read-only restriction.

    :param h5file: The open HDF5 file handle.
    :type h5file: h5py.File
    :param collection_spec: The collection specification defining available datasets.
    :type collection_spec: DataCollectionSpec
    :param sim: The simulation index.
    :type sim: int
    """

    def __init__(self, h5file, collection_spec, sim):
        self._h5file = h5file
        self._spec = collection_spec
        self._sim = sim

    @property
    def datasets(self):
        """List of dataset names available for this simulation.

        :return: Dataset names from the collection spec (excluding those with
                 ``data_location=None``).
        :rtype: list[str]
        """
        return [
            name
            for name, ds_spec in self._spec.data.items()
            if ds_spec.data_location is not None
        ]

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._spec.data:
            ds_spec = self._spec.data[name]
            if ds_spec.data_location is None:
                raise AttributeError(
                    f"Dataset '{name}' is a derived dataset and not stored on disk."
                )
            path = ds_spec.data_location.format(sim=self._sim)
            return self._h5file[path]
        raise AttributeError(
            f"'{type(self).__name__}' has no dataset '{name}'. "
            f"Available datasets: {self.datasets}"
        )

    def __repr__(self):
        return (
            f"<SimulationView sim={self._sim}, "
            f"datasets={len(self.datasets)}>"
        )


class DataCollection:
    """Structured interface to one data collection within the HDF5 file.

    Whether datasets are writable depends on the mode of the underlying
    :class:`h5py.File`; this class imposes no read-only restriction.

    Behavior depends on whether simulations are combined:

    **Combined** (``simulations_combined=True``, e.g., microscale_out):
    Dot-access to dataset names returns the ``h5py.Dataset`` object
    (lazy-loaded — data is not read until sliced).

    **Per-simulation** (``simulations_combined=False``, e.g., macroscale_out):
    Use ``collection[sim_index]`` to get a :class:`SimulationView`, then
    access datasets on it. Direct dot-access to dataset names raises
    ``TypeError`` with a helpful message.

    :param name: The collection name (e.g., ``"microscale_out"``).
    :type name: str
    :param h5file: The open HDF5 file handle.
    :type h5file: h5py.File
    :param collection_spec: The spec defining this collection's datasets.
    :type collection_spec: DataCollectionSpec
    :param num_sims: Number of simulations. Required for per-simulation
        collections (``simulations_combined=False``). Sourced from
        ``MacroParameters.macro_simulations``.
    :type num_sims: int | None
    """

    def __init__(self, name, h5file, collection_spec, num_sims=None):
        self._name = name
        self._h5file = h5file
        self._spec = collection_spec
        if not collection_spec.simulations_combined:
            if num_sims is None:
                raise ValueError(
                    f"Per-simulation collection '{name}' requires num_sims."
                )
            self._num_sims = num_sims

    @property
    def datasets(self):
        """List of dataset names in this collection.

        :return: Dataset names from the spec (excluding those with
                 ``data_location=None``).
        :rtype: list[str]
        """
        return [
            name
            for name, ds_spec in self._spec.data.items()
            if ds_spec.data_location is not None
        ]

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._spec.data:
            if not self._spec.simulations_combined:
                raise TypeError(
                    f"Collection '{self._name}' stores data per-simulation. "
                    f"Index by simulation first: "
                    f"{self._name}[sim_index].{name}"
                )
            ds_spec = self._spec.data[name]
            if ds_spec.data_location is None:
                raise AttributeError(
                    f"Dataset '{name}' is a derived dataset and not stored on disk."
                )
            return self._h5file[ds_spec.data_location]
        raise AttributeError(
            f"'{type(self).__name__}' has no dataset '{name}'. "
            f"Available datasets: {self.datasets}"
        )

    def __getitem__(self, index):
        if self._spec.simulations_combined:
            raise TypeError(
                f"Collection '{self._name}' has combined simulations "
                f"and does not support indexing. "
                f"Access datasets directly: {self._name}.dataset_name"
            )
        if not isinstance(index, int):
            raise TypeError(
                f"Simulation index must be an integer, got {type(index).__name__}"
            )
        if index < 0 or index >= self._num_sims:
            raise IndexError(
                f"Simulation index {index} out of range "
                f"(0 to {self._num_sims - 1})"
            )
        return SimulationView(self._h5file, self._spec, index)

    def __len__(self):
        if self._spec.simulations_combined:
            raise TypeError(
                f"Collection '{self._name}' has combined simulations "
                f"and does not have a length. "
                f"Access datasets directly: {self._name}.dataset_name"
            )
        return self._num_sims

    def __contains__(self, name):
        return name in self._spec.data and self._spec.data[name].data_location is not None

    def __repr__(self):
        kind = "combined" if self._spec.simulations_combined else "per-simulation"
        parts = f"<DataCollection '{self._name}', {kind}, datasets={len(self.datasets)}"
        if not self._spec.simulations_combined:
            parts += f", simulations={self._num_sims}"
        parts += ">"
        return parts


class DerivedDataCollection:
    """Read-only interface to a lazily-generated derived data collection.

    Provides the same dot-access interface as :class:`DataCollection` for
    combined collections, but stores numpy arrays in memory rather than
    ``h5py.Dataset`` objects.  Data is generated on first access by calling
    the supplied *generator* callable, cached for subsequent accesses, and
    immediately marked non-writeable — attempts to modify a returned array
    will raise :exc:`ValueError`.

    :param name: The collection name (e.g., ``"macroscale_in"``).
    :type name: str
    :param collection_spec: The spec defining this collection's datasets.
    :type collection_spec: DataCollectionSpec
    :param generator: A callable that returns a :class:`DataCollectionType`
        dict mapping dataset names (and ``"params"``) to numpy arrays.
    :type generator: callable
    """

    def __init__(self, name, collection_spec, generator):
        self._name = name
        self._spec = collection_spec
        self._generator = generator
        self._data = None  # populated on first access

    def _ensure_generated(self):
        """Call the generator if data has not yet been generated.

        After generation, all numpy arrays in the result are marked
        non-writeable so callers cannot mutate the cached data.
        """
        if self._data is None:
            self._data = self._generator()
            for value in self._data.values():
                if isinstance(value, np.ndarray):
                    value.flags.writeable = False

    @property
    def datasets(self):
        """List of dataset names in this collection.

        Unlike :attr:`DataCollection.datasets`, this includes *all* datasets
        from the spec regardless of ``data_location``, since derived
        collections have no on-disk storage.

        :return: All dataset names from the spec.
        :rtype: list[str]
        """
        return list(self._spec.data.keys())

    @property
    def params(self):
        """Parameters dict computed during generation.

        Triggers generation if data has not yet been generated.

        :return: The ``"params"`` entry from the generator output.
        :rtype: dict
        """
        self._ensure_generated()
        return self._data.get("params")

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._spec.data:
            self._ensure_generated()
            if name in self._data:
                return self._data[name]
            raise AttributeError(
                f"Dataset '{name}' is in the spec but was not produced "
                f"by the generator."
            )
        raise AttributeError(
            f"'{type(self).__name__}' has no dataset '{name}'. "
            f"Available datasets: {self.datasets}"
        )

    def __contains__(self, name):
        return name in self._spec.data

    def __repr__(self):
        status = "generated" if self._data is not None else "pending"
        return (
            f"<DerivedDataCollection '{self._name}', "
            f"status={status}, datasets={len(self.datasets)}>"
        )


class DataStore:
    """Interface to an HDF5 file following the v2.0.0 data specification.

    Provides dot-access to data collections and parameters::

        with DataStore("run01", "/path/to/data") as ds:
            ds.microscale_out.tpa_leaving_time[:10]  # lazy h5py.Dataset
            ds.macroscale_out[3].fiber_degrade_time   # per-simulation access
            ds.micro_params.fiber_radius              # loaded MicroParameters
            ds.macro_params.rows                      # loaded MacroParameters

    To create a new DataStore, use the :meth:`create` class method::

        ds = DataStore.create("run01", "/path/to/data", micro_params)

    To add macroscale data to an existing microscale DataStore
    (modifies the DataStore in place; requires writable mode)::

        ds.initialize_macroscale(macro_params)

    :param run_code: The run identifier (used to construct the HDF5 filename).
    :type run_code: str
    :param path: Directory containing the HDF5 file.
    :type path: str
    :param mode: HDF5 file mode passed directly to :class:`h5py.File`.
        ``"r"`` (default) opens read-only (file must exist).
        ``"a"`` opens read/write (creates file if it doesn't exist).
        ``"w"`` creates a new file (truncates if it exists).
    :type mode: str
    """

    def __init__(self, run_code, path, mode="r"):
        self._run_code = run_code
        self._path = path
        self._mode = mode
        self._hdf5_path = os.path.join(path, f"{run_code}.h5")
        self._file = h5py.File(self._hdf5_path, mode)

        # Validate dataspec version (also warns if CONVERTED_FROM_ATTR is set)
        try:
            _validate_hdf5_version(
                self._file, self._hdf5_path, COMPATIBLE_DATASPEC_VERSION
            )
        except ValueError:
            self.close()
            raise
        self._dataspec_version = COMPATIBLE_DATASPEC_VERSION

        self._collections = {}
        self._micro_params = None
        self._macro_params = None

        spec = dataspec[COMPATIBLE_DATASPEC_VERSION]

        # Detect which HDF5 groups are present
        present_groups = {}
        for coll_name, coll_spec in spec.items():
            group_path = _group_path_from_spec(coll_spec)
            if group_path is None:
                # Derived collection (e.g., macroscale_in) — skip
                continue
            if group_path in self._file:
                present_groups[coll_name] = coll_spec

        # Validate collection dependencies
        if "macroscale_out" in present_groups and "microscale_out" not in present_groups:
            self.close()
            raise ValueError(
                "HDF5 file contains macroscale_out data without microscale_out data. "
                "Macroscale output requires microscale output to be present."
            )

        # Load parameters using fileops.read_dataset()
        # Parameters are required whenever the corresponding data group exists.
        if "microscale_out" in present_groups:
            micro_spec = spec["microscale_out"]
            raw_params = read_dataset(self._hdf5_path, micro_spec.params)
            inner = raw_params.get(micro_spec.params.data_location, {})
            if not inner:
                self.close()
                raise ValueError(
                    "microscale_out data is present but the micro_data group "
                    "has no parameter attributes."
                )
            self._micro_params = MicroParameters.parse_from_basedict(inner)

        if "macroscale_out" in present_groups:
            macro_spec = spec["macroscale_out"]
            raw_params = read_dataset(self._hdf5_path, macro_spec.params)
            inner = raw_params.get(macro_spec.params.data_location, {})
            if not inner:
                self.close()
                raise ValueError(
                    "macroscale_out data is present but the macro_data group "
                    "has no parameter attributes."
                )
            inner["micro_params"] = self._micro_params
            self._macro_params = MacroParameters.parse_from_basedict(inner)

        # Create DataCollection objects (after params are loaded)
        for coll_name, coll_spec in present_groups.items():
            num_sims = None
            if not coll_spec.simulations_combined:
                num_sims = self._macro_params.macro_simulations
            self._collections[coll_name] = DataCollection(
                coll_name, self._file, coll_spec, num_sims=num_sims
            )

        # Register derived collections (lazy generation)
        if self._can_generate_macroscale_in():
            macro_in_spec = spec["macroscale_in"]
            self._collections["macroscale_in"] = DerivedDataCollection(
                "macroscale_in",
                macro_in_spec,
                self._build_macroscale_in_generator(),
            )

        self._status = DataStatus.INITIALIZED

    @property
    def mode(self):
        """The file mode this DataStore was opened with.

        :rtype: str
        """
        return self._mode

    @property
    def dataspec_version(self):
        """The dataspec version string read from the HDF5 file.

        :rtype: str
        """
        return self._dataspec_version

    @property
    def collections(self):
        """Dictionary of available data collections.

        :return: Mapping of collection name to :class:`DataCollection`.
        :rtype: dict[str, DataCollection]
        """
        return dict(self._collections)

    @property
    def micro_params(self):
        """Microscale parameters loaded from HDF5 attributes.

        ``None`` only when no microscale data is present.

        :rtype: MicroParameters | None
        """
        return self._micro_params

    @property
    def macro_params(self):
        """Macroscale parameters loaded from HDF5 attributes.

        ``None`` only when no macroscale data is present.

        :rtype: MacroParameters | None
        """
        return self._macro_params

    @property
    def status(self):
        """The current state of this DataStore.

        :rtype: DataStatus
        """
        return self._status

    @property
    def hdf5_state(self) -> HDF5State:
        """Pipeline lifecycle state inferred from dataset fill status.

        Inspects the ``micro_data/tpa_leaving_time`` dataset to determine
        whether microscale results have been written, and (when the macroscale
        group is present) ``macro_data/sim_00/snapshot_time`` to determine
        whether macroscale results have been written.

        :return: One of :class:`HDF5State` ``MICRO_EMPTY``, ``MICRO_FILLED``,
            ``MACRO_EMPTY``, ``MACRO_FILLED``, or ``INCONSISTENT``.
        :rtype: HDF5State
        """
        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
        check_loc = micro_spec.data["tpa_leaving_time"].data_location
        micro_filled = (
            check_loc in self._file and self._file[check_loc].shape[0] > 0
        )

        if "macroscale_out" not in self._collections:
            return HDF5State.MICRO_FILLED if micro_filled else HDF5State.MICRO_EMPTY

        if not micro_filled:
            return HDF5State.INCONSISTENT

        macro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["macroscale_out"]
        check_loc = macro_spec.data["snapshot_time"].data_location.format(sim=0)
        macro_filled = (
            check_loc in self._file and self._file[check_loc].shape[0] > 0
        )
        return HDF5State.MACRO_FILLED if macro_filled else HDF5State.MACRO_EMPTY

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in ("mode", "micro_params", "macro_params", "collections", "status", "hdf5_state"):
            # These are properties — if we're here, the object isn't
            # fully initialized yet. Avoid infinite recursion.
            raise AttributeError(name)
        # Check collections
        try:
            collections = object.__getattribute__(self, "_collections")
        except AttributeError:
            raise AttributeError(name)
        if name in collections:
            return collections[name]
        raise AttributeError(
            f"'{type(self).__name__}' has no attribute '{name}'. "
            f"Available collections: {list(collections.keys())}"
        )

    def __str__(self):
        """Print the datastore HDF5 structure."""
        return h5_tree(self._file)

    def __repr__(self):
        colls = list(self._collections.keys())
        return f"<DataStore '{self._run_code}', collections={colls}>"

    @staticmethod
    def _create_empty_datasets(hdf5_path, collection_spec, params=None,
                               num_sims=None):
        """Create zero-length HDF5 datasets for all storable datasets in a collection.

        Each dataset is created with shape zero in variable dimensions (``-1``
        in the spec) and ``maxshape=None`` for those dimensions so they can be
        resized later when data is written.

        :param hdf5_path: Path to the HDF5 file (must already exist).
        :type hdf5_path: str
        :param collection_spec: The collection specification.
        :type collection_spec: DataCollectionSpec
        :param params: Parameter dictionary for resolving dynamic shape
            dimensions, in ``{"micro_params": {...}, "macro_params": {...}}``
            format. Only needed if spec shapes contain string references.
        :type params: dict, optional
        :param num_sims: Number of simulations. Required for per-simulation
            collections (``simulations_combined=False``).
        :type num_sims: int, optional
        """
        with h5py.File(hdf5_path, "a") as f:
            for ds_name, ds_spec in collection_spec.data.items():
                if ds_spec.data_location is None:
                    # Derived dataset — not stored on disk
                    continue

                # Resolve the spec shape (replace string refs with param values)
                resolved = parse_shape(ds_spec.shape, params)

                # Build initial (empty) shape and maxshape
                empty_shape = tuple(
                    0 if dim == -1 else dim for dim in resolved
                )
                maxshape = tuple(
                    None if dim == -1 else dim for dim in resolved
                )

                if collection_spec.simulations_combined:
                    # Combined collection: one dataset at ds_spec.data_location
                    f.create_dataset(
                        ds_spec.data_location,
                        shape=empty_shape,
                        maxshape=maxshape,
                        dtype=ds_spec.dtype,
                        compression="gzip",
                    )
                else:
                    # Per-simulation: create dataset in each sim group
                    for sim in range(num_sims):
                        path = ds_spec.data_location.format(sim=sim)
                        f.create_dataset(
                            path,
                            shape=empty_shape,
                            maxshape=maxshape,
                            dtype=ds_spec.dtype,
                            compression="gzip",
                        )

    @classmethod
    def create(cls, run_code, path, micro_params, force=False):
        """Create a new DataStore with microscale parameters and empty datasets.

        Creates a new HDF5 file, writes the dataspec version attribute,
        stores ``micro_params`` as HDF5 attributes on the ``micro_data``
        group, and creates zero-length datasets for all microscale_out
        datasets defined in the v2.0.0 specification.

        Microscale ``init_*`` provenance is stamped on the ``micro_data``
        params group before returning (via :meth:`stamp_provenance`), so the
        returned file already presents a clean "initialized" state.

        :param run_code: The Run identifier (used to construct the HDF5
            filename as ``{run_code}.h5``).
        :type run_code: str
        :param path: Directory where the HDF5 file will be created.
        :type path: str
        :param micro_params: The microscale parameters to store.
        :type micro_params: MicroParameters
        :param force: If ``True`` and the target HDF5 file already exists,
            delete it and recreate from scratch.  Any previous execution
            provenance attributes on the file are lost along with the file.
        :type force: bool
        :return: A new DataStore opened in ``"a"`` (read/write) mode with
            microscale_out collection.
        :rtype: DataStore
        :raises FileExistsError: If the HDF5 file already exists and
            ``force`` is ``False``.
        :raises TypeError: If ``micro_params`` is not a
            :class:`~lysis.config.parameters.MicroParameters` instance.
        """
        if not isinstance(micro_params, MicroParameters):
            raise TypeError(
                f"Expected MicroParameters, got {type(micro_params).__name__}"
            )

        hdf5_path = os.path.join(path, f"{run_code}.h5")
        if os.path.exists(hdf5_path):
            if not force:
                raise FileExistsError(f"HDF5 file already exists: {hdf5_path}")
            os.remove(hdf5_path)

        spec = dataspec[COMPATIBLE_DATASPEC_VERSION]
        micro_spec = spec["microscale_out"]

        # Write micro parameters as HDF5 attributes on the micro_data group.
        # write_dataset -> _write_hdf5_attr -> ensure_hdf5_version (creates file)
        param_group_name = micro_spec.params.data_location.replace(
            "_data", "_params"
        )
        params_data = {param_group_name: micro_params.to_basedict()}
        write_dataset(params_data, hdf5_path, micro_spec.params)

        # Create empty datasets for microscale_out
        params_dict = {"micro_params": micro_params.to_basedict()}
        cls._create_empty_datasets(hdf5_path, micro_spec, params=params_dict)

        # Stamp microscale init provenance so the freshly-created file already
        # presents a clean "initialized" state.  Both the init-experiment CLI
        # and the import-failure rollback reach the stamp through this path.
        ds = cls(run_code, path, mode="a")
        ds.stamp_provenance("micro", "init")
        return ds

    @classmethod
    def rename(cls, old_run_code, new_run_code, path):
        """Rename the HDF5 file for a Run and record the previous name.

        Renames ``{path}/{old_run_code}.h5`` to
        ``{path}/{new_run_code}.h5`` on disk, then writes the previous
        ``run_code`` into the HDF5 root attribute
        ``CONST.RENAMED_FROM_ATTR``.  If that attribute already contains a
        chronological history (earlier names joined by ``" -> "``), the old
        ``run_code`` is appended to the end so the attribute always reads
        oldest → most recent previous name.

        :param old_run_code: Current ``run_code`` — stem of the existing
            HDF5 filename.
        :type old_run_code: str
        :param new_run_code: New ``run_code`` — stem of the destination
            HDF5 filename.
        :type new_run_code: str
        :param path: Directory containing the HDF5 file.
        :type path: str
        :raises ValueError: If ``new_run_code`` equals ``old_run_code``.
        :raises FileNotFoundError: If the source HDF5 file does not exist.
        :raises FileExistsError: If an HDF5 file already exists at the
            destination path.
        """
        if new_run_code == old_run_code:
            raise ValueError(
                f"new_run_code is identical to old_run_code: {new_run_code!r}"
            )

        old_path = os.path.join(path, f"{old_run_code}.h5")
        new_path = os.path.join(path, f"{new_run_code}.h5")

        if not os.path.isfile(old_path):
            raise FileNotFoundError(f"HDF5 file not found: {old_path}")
        if os.path.exists(new_path):
            raise FileExistsError(f"HDF5 file already exists: {new_path}")

        os.rename(old_path, new_path)

        attr = CONST.RENAMED_FROM_ATTR
        with h5py.File(new_path, "a") as f:
            existing = f.attrs.get(attr)
            if existing is None:
                history = old_run_code
            else:
                history = f"{existing} -> {old_run_code}"
            f.attrs[attr] = history

    def initialize_macroscale(self, macro_params, force=False):
        """Add macroscale parameters and empty datasets to this DataStore.

        Reads the completed microscale output to compute
        :attr:`~lysis.config.parameters.MacroParameters.forced_unbind`
        automatically from the ``tpa_unbound_by_pli`` and
        ``tpa_unbound_kinetic`` datasets, then writes the full macroscale
        structure (parameters + empty per-simulation dataset stubs) to the
        HDF5 file.

        **This method must be called after all microscale Simulations have
        completed and their results have been written to the DataStore.**
        Calling it before microscale output exists will raise
        :exc:`ValueError`.

        The value of ``macro_params.forced_unbind`` is silently replaced
        with the value computed from microscale data; do not rely on
        whatever value was passed in.

        Macroscale ``init_*`` provenance is stamped on the ``macro_data``
        params group once the empty structure is in place (via
        :meth:`stamp_provenance`), so the file presents a clean macroscale
        "initialized" state.  Both the init-macroscale CLI and the
        import-failure rollback reach the stamp through this path.

        Modifies the DataStore **in place** and returns ``None``.
        Requires the DataStore to be opened in a writable mode
        (e.g., ``"a"``)::

            with DataStore(run_code, path, mode="a") as ds:
                ds.initialize_macroscale(run.macro_params)

        :param macro_params: The macroscale parameters to store.  All
            fields except ``forced_unbind`` are used as provided.
        :type macro_params: MacroParameters
        :param force: If ``True`` and ``macroscale_out`` already exists,
            wipe the ``macro_data`` group and any ``log_files/macro_log__sim_*``
            datasets before re-initialising.  Any previously stored
            execution provenance attributes are lost along with the group.
        :type force: bool
        :raises IOError: If the DataStore is opened in read-only mode.
        :raises TypeError: If ``macro_params`` is not a
            :class:`~lysis.config.parameters.MacroParameters` instance.
        :raises ValueError: If microscale_out is not present, if
            macroscale_out is already present and ``force`` is ``False``,
            or if the microscale datasets are empty (i.e. no Simulations
            have been written).
        """
        if self._mode == "r":
            raise IOError(
                "Cannot initialize macroscale on a read-only DataStore. "
                "Open with mode='a' or use DataStore.create()."
            )
        if not isinstance(macro_params, MacroParameters):
            raise TypeError(
                f"Expected MacroParameters, got {type(macro_params).__name__}"
            )
        if "microscale_out" not in self._collections:
            raise ValueError(
                "Cannot add macroscale data: microscale_out is not present. "
                "Use DataStore.create() first."
            )
        if "macroscale_out" in self._collections:
            if not force:
                raise ValueError(
                    "macroscale_out is already present in this DataStore."
                )
            # Wipe the existing macro_data group and per-sim macro log
            # datasets. Execution provenance attrs live on macro_data and
            # are cleared automatically when the group is deleted.
            if "macro_data" in self._file:
                del self._file["macro_data"]
            if "log_files" in self._file:
                macro_log_keys = [
                    k for k in self._file["log_files"]
                    if k.startswith("macro_log__sim_")
                ]
                for key in macro_log_keys:
                    del self._file[f"log_files/{key}"]
            self._file.flush()

        # Read microscale unbinding arrays and compute forced_unbind from data.
        # This must happen before self._file is closed.
        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
        pli_loc = micro_spec.data["tpa_unbound_by_pli"].data_location
        kin_loc = micro_spec.data["tpa_unbound_kinetic"].data_location
        tpa_unbound_by_pli = self._file[pli_loc][:]
        tpa_unbound_kinetic = self._file[kin_loc][:]

        if tpa_unbound_by_pli.shape[0] == 0:
            raise ValueError(
                "Cannot initialize macroscale: microscale datasets are empty. "
                "Ensure all microscale Simulations have completed before "
                "calling initialize_macroscale()."
            )

        computed_forced_unbind = MacroParameters.calculate_forced_unbind(
            tpa_unbound_by_pli, tpa_unbound_kinetic
        )
        macro_params = dataclasses.replace(
            macro_params, forced_unbind=computed_forced_unbind
        )

        spec = dataspec[COMPATIBLE_DATASPEC_VERSION]
        macro_spec = spec["macroscale_out"]
        hdf5_path = self._hdf5_path
        micro_params = self._micro_params

        # Close the file handle before writing (fileops opens its own handles)
        self._file.close()

        # Write macro parameters as HDF5 attributes on the macro_data group
        param_group_name = macro_spec.params.data_location.replace(
            "_data", "_params"
        )
        params_data = {param_group_name: macro_params.to_basedict()}
        write_dataset(params_data, hdf5_path, macro_spec.params)

        # Create empty per-simulation datasets for macroscale_out
        params_dict = {
            "micro_params": micro_params.to_basedict(),
            "macro_params": macro_params.to_basedict(),
        }
        num_sims = macro_params.macro_simulations
        type(self)._create_empty_datasets(
            hdf5_path, macro_spec, params=params_dict, num_sims=num_sims
        )

        # Re-initialize in place (reloads all collections, params, etc.)
        self.__init__(self._run_code, self._path, mode=self._mode)

        # Stamp macroscale init provenance now that the empty structure is in
        # place, mirroring create()'s micro stamp.
        self.stamp_provenance("macro", "init")

    # ------------------------------------------------------------------
    #  Provenance stamping
    # ------------------------------------------------------------------

    def _scale_to_params_location(self, scale):
        """Resolve a ``"micro"``/``"macro"`` scale name to the HDF5 params
        group location, via the v2.0.0 dataspec.

        :raises ValueError: If *scale* is not ``"micro"`` or ``"macro"``.
        """
        if scale == "micro":
            collection_name = "microscale_out"
        elif scale == "macro":
            collection_name = "macroscale_out"
        else:
            raise ValueError(
                f"scale must be 'micro' or 'macro'; got {scale!r}"
            )
        spec = dataspec[COMPATIBLE_DATASPEC_VERSION][collection_name]
        return spec.params.data_location

    def stamp_provenance(
        self,
        scale,
        kind,
        *,
        executable=None,
        backend_override=None,
        replace_backend_attrs=None,
    ):
        """Stamp provenance attributes onto the per-scale params group.

        Writes the ``init_*``, ``pipeline_*``, or ``backend_*`` attribute
        family to the HDF5 group at ``{scale}_data`` (the params group
        for the named scale).

        :param scale: ``"micro"`` (targets ``micro_data``) or ``"macro"``
            (targets ``macro_data``).
        :type scale: str
        :param kind: One of ``"init"``, ``"pipeline"``, ``"backend"``.

            * ``"init"`` calls
              :func:`~lysis.tools.provenance.gather_init_provenance`.
            * ``"pipeline"`` calls
              :func:`~lysis.tools.provenance.gather_pipeline_provenance`.
            * ``"backend"`` calls
              :func:`~lysis.tools.provenance.gather_backend_provenance`
              against *executable*, stamps ``backend_type = "fortran"``,
              and additionally writes any keys in *backend_override*.
        :type kind: str
        :param executable: Required when ``kind="backend"`` — path to the
            Fortran binary whose ``--version`` output supplies commit /
            dirty / compiler stamps.  May be ``None`` when
            *replace_backend_attrs* is supplied (historical-build mode,
            where the binary's commit/dirty/compiler are synthesised
            externally rather than queried).  Ignored for other kinds.
        :type executable: pathlib.Path or str, optional
        :param backend_override: Optional dict of extra attrs to merge in
            *on top of* the gathered (or replaced) attrs — typically
            ``{stale_backend_override: True}`` from
            :func:`~lysis.tools.provenance.verify_binary_matches_source`
            when the staleness check was overridden.  Ignored when
            ``kind != "backend"``.
        :type backend_override: dict, optional
        :param replace_backend_attrs: Optional pre-computed dict that
            *replaces* the default
            :func:`~lysis.tools.provenance.gather_backend_provenance`
            call when ``kind="backend"``.  Used by the historical-build
            workflow to ship a synthesised provenance dict (commit SHA,
            ``iso_fortran_env`` compiler probe output, and the
            ``backend_historical = True`` marker) for a binary whose own
            ``--version`` may not match — or be absent.  When supplied,
            *executable* may be ``None`` and no subprocess is spawned
            against the binary.  ``backend_override`` (if also supplied)
            is merged on top.  Ignored for other kinds.
        :type replace_backend_attrs: dict, optional
        :raises IOError: If the DataStore is in read-only mode.
        :raises ValueError: For an unknown *scale*/*kind*, a missing
            target group, or ``kind="backend"`` with neither
            *executable* nor *replace_backend_attrs*.
        """
        if self._mode == "r":
            raise IOError(
                "Cannot stamp provenance on a read-only DataStore. "
                "Open with mode='a'."
            )

        params_location = self._scale_to_params_location(scale)
        if params_location not in self._file:
            raise ValueError(
                f"Cannot stamp provenance: params group "
                f"'{params_location}' does not exist in this DataStore.  "
                f"Initialise the {scale}scale collection first."
            )

        if kind == "init":
            from ..tools.provenance import gather_init_provenance  # noqa: PLC0415
            attrs = gather_init_provenance()
        elif kind == "pipeline":
            from ..tools.provenance import gather_pipeline_provenance  # noqa: PLC0415
            attrs = gather_pipeline_provenance()
        elif kind == "backend":
            if replace_backend_attrs is not None:
                attrs = dict(replace_backend_attrs)
            else:
                if executable is None:
                    raise ValueError(
                        "stamp_provenance(kind='backend') requires either "
                        "an 'executable' argument or 'replace_backend_attrs'."
                    )
                from ..tools.provenance import gather_backend_provenance  # noqa: PLC0415
                attrs = gather_backend_provenance(executable)
            # Every current run is Fortran; the future Python backend will
            # stamp "python" here once wired into run-* (see #35).
            attrs = {**attrs, CONST.BACKEND_TYPE_ATTR: "fortran"}
            if backend_override:
                attrs = {**attrs, **backend_override}
        else:
            raise ValueError(
                f"kind must be 'init', 'pipeline', or 'backend'; "
                f"got {kind!r}"
            )

        group = self._file[params_location]
        for attr_name, attr_value in attrs.items():
            group.attrs[attr_name] = attr_value

    def read_init_provenance(self, scale):
        """Read the ``init_*`` attributes stamped on the per-scale params group.

        :param scale: ``"micro"`` or ``"macro"``.
        :type scale: str
        :return: Dict with keys :data:`CONST.INIT_VERSION_ATTR`,
            :data:`CONST.INIT_DIRTY_ATTR`,
            :data:`CONST.INIT_TIMESTAMP_ATTR`,
            :data:`CONST.INIT_HOSTNAME_ATTR` — or ``None`` if the params
            group exists but the file pre-dates the init-stamp feature
            (no ``init_version`` attribute present).
        :rtype: dict or None
        :raises ValueError: For an unknown *scale* or a missing params
            group.
        """
        params_location = self._scale_to_params_location(scale)
        if params_location not in self._file:
            raise ValueError(
                f"Cannot read init provenance: params group "
                f"'{params_location}' does not exist in this DataStore."
            )
        attrs = self._file[params_location].attrs
        if CONST.INIT_VERSION_ATTR not in attrs:
            return None
        keys = (
            CONST.INIT_VERSION_ATTR,
            CONST.INIT_DIRTY_ATTR,
            CONST.INIT_TIMESTAMP_ATTR,
            CONST.INIT_HOSTNAME_ATTR,
        )
        result = {}
        for k in keys:
            if k in attrs:
                val = attrs[k]
                if isinstance(val, bytes):
                    val = val.decode("utf-8")
                result[k] = val
        return result

    def _write_array_to_dataset(self, ds_spec, arr, sim=None):
        """Write a numpy array to an existing empty HDF5 dataset.

        Handles both numeric arrays (resize then slice-assign) and variable-
        length string arrays (delete the empty placeholder and recreate with
        data, since h5py string datasets cannot be reliably resized and
        slice-assigned).

        :param ds_spec: Specification for the target dataset.
        :type ds_spec: DataSetSpec
        :param arr: The array to write.
        :type arr: numpy.ndarray
        :param sim: Simulation index for per-simulation datasets.  ``None``
            for combined-simulation datasets.
        :type sim: int, optional
        """
        if sim is not None:
            loc = ds_spec.data_location.format(sim=sim)
        else:
            loc = ds_spec.data_location

        if ds_spec.dtype == h5py.string_dtype():
            # Variable-length string datasets: delete the shape-0 placeholder
            # and recreate with the actual data.
            if loc in self._file:
                del self._file[loc]
            self._file.create_dataset(
                loc,
                data=arr,
                dtype=h5py.string_dtype(),
                maxshape=(None,),
                compression="gzip",
            )
        else:
            h5ds = self._file[loc]
            h5ds.resize(arr.shape)
            h5ds[:] = arr

    def import_collection(
        self,
        collection_name,
        source_spec,
        source_path,
        file_codes,
        param_overrides=None,
        param_aliases=None,
        backend_executable=None,
        backend_override=None,
        replace_backend_attrs=None,
    ):
        """Import a data collection from an external source into this DataStore.

        Reads data from *source_path* using the *source_spec* format, converts
        it to v2.0.0 if necessary, then writes the resulting numpy arrays
        directly into the existing empty HDF5 datasets managed by this
        DataStore.

        The target collection must already exist with empty datasets:

        - For ``"microscale_out"``: call :meth:`create` first.
        - For ``"macroscale_out"``: call :meth:`initialize_macroscale` first
          (which itself requires filled microscale data).

        Modifies the DataStore **in place** and returns ``None``.
        Requires the DataStore to be opened in writable mode (``"a"``)::

            ds = DataStore.create(run_code, path, micro_params)
            ds.import_collection(
                "microscale_out", "v1.95.0", fortran_path,
                [micro_file_code], param_overrides=overrides
            )

        :param collection_name: The collection to populate.  Must be
            ``"microscale_out"`` or ``"macroscale_out"``.
        :type collection_name: str
        :param source_spec: Dataspec version of the source data, or a tag
            alias (e.g. ``"v1.95.0"``, ``"fortran"``, ``"v2.0.0"``,
            ``"hdf5"``).  Tag aliases are resolved automatically.
        :type source_spec: str
        :param source_path: Path to the source data.  For Fortran-format
            specs, a directory containing simulation output files.  For
            HDF5-format specs, the full path to the HDF5 file.
        :type source_path: str
        :param file_codes: List of file-code strings passed to
            :func:`~lysis.dataio.fileops.read_data_collection`.  Typically
            one entry (e.g. ``["_PLG2_tPA01_TB-xiii"]`` for a Fortran
            microscale source, or ``[""]`` for an HDF5 source).
        :type file_codes: list[str]
        :param param_overrides: Optional key-value pairs injected into the
            read parameters before conversion.  Useful when Fortran log files
            omit fields that downstream converters require.
        :type param_overrides: dict, optional
        :param param_aliases: Optional ``{python_name: fortran_name}`` mapping
            used to rename Fortran parameter keys before conversion.  For
            example, ``{"micro_simulations": "runs"}`` renames the ``runs``
            key produced by v1.90.0 log parsing.
        :type param_aliases: dict, optional
        :param backend_executable: Path to the Fortran binary that produced
            the source data.  When provided, the backend's commit/dirty/
            compiler stamps (from ``<executable> --version``) are written
            to the per-scale params group via
            :meth:`stamp_provenance`.  ``None`` (the default) skips the
            backend stamp — appropriate for HDF5→HDF5 conversions and tests.
        :type backend_executable: pathlib.Path or str, optional
        :param backend_override: Optional dict of override-only attrs from
            :func:`~lysis.tools.provenance.verify_binary_matches_source`
            (typically ``{stale_backend_override: True}``) when the binary
            preflight check was bypassed.  Merged into the backend-stamp
            group attrs.  Ignored when neither ``backend_executable`` nor
            ``replace_backend_attrs`` is provided.
        :type backend_override: dict, optional
        :param replace_backend_attrs: Optional pre-computed backend-provenance
            dict from
            :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
            (historical-build workflow).  When supplied, the backend stamp
            is written from this dict instead of querying
            ``<backend_executable> --version``; ``backend_executable`` may
            be ``None`` in that case.  See
            :meth:`stamp_provenance` for details.
        :type replace_backend_attrs: dict, optional
        :raises IOError: If the DataStore is in read-only mode.
        :raises ValueError: If *collection_name* is not ``"microscale_out"``
            or ``"macroscale_out"``, if the target collection does not yet
            exist (call :meth:`create` / :meth:`initialize_macroscale` first),
            if any target dataset already contains data, or if the number of
            simulations in the source does not match
            ``macro_params.macro_simulations`` (for ``"macroscale_out"``).
        :raises ImportCollectionError: If the import fails *after* HDF5
            mutation has begun (e.g. a dataset write errors).  Before
            re-raising, the target collection is rolled back to its empty,
            freshly-initialized state (``MICRO_EMPTY`` / ``MACRO_EMPTY``) so
            the file never presents half-ingested data as filled; the original
            failure is chained as ``__cause__``.

        .. note::

           On failure during the write phase the target collection is reverted
           to empty via the same machinery as the ``init-*`` pathway, so
           downstream code never processes a partially-imported collection.
           Any on-disk source (Fortran output / logs) is left uningested.
        """
        # ------------------------------------------------------------------
        # Precondition checks
        # ------------------------------------------------------------------
        if self._mode == "r":
            raise IOError(
                "Cannot import data on a read-only DataStore. "
                "Open with mode='a'."
            )
        if collection_name not in ("microscale_out", "macroscale_out"):
            raise ValueError(
                f"Cannot import collection '{collection_name}'. "
                "Only 'microscale_out' and 'macroscale_out' are importable."
            )
        if collection_name not in self._collections:
            if collection_name == "microscale_out":
                raise ValueError(
                    "Cannot import microscale_out: collection does not exist. "
                    "Call DataStore.create() first."
                )
            else:
                raise ValueError(
                    "Cannot import macroscale_out: collection does not exist. "
                    "Call initialize_macroscale() first."
                )

        # All storable datasets in the target collection must be empty.
        target_spec = dataspec[COMPATIBLE_DATASPEC_VERSION][collection_name]
        for ds_name, ds_spec in target_spec.data.items():
            if ds_spec.data_location is None:
                continue
            if target_spec.simulations_combined:
                loc = ds_spec.data_location
            else:
                loc = ds_spec.data_location.format(sim=0)
            if loc in self._file and self._file[loc].shape[0] != 0:
                raise ValueError(
                    f"Cannot import {collection_name}: dataset '{ds_name}' "
                    "already contains data. Import is only allowed once, "
                    "on empty datasets."
                )

        # ------------------------------------------------------------------
        # Resolve tag aliases in source_spec
        # ------------------------------------------------------------------
        source_version = source_spec
        while source_version in tags:
            source_version = tags[source_version]

        # ------------------------------------------------------------------
        # Read from the external source
        # ------------------------------------------------------------------
        # Import convert_data lazily to avoid circular imports
        # (dataconvert → geometry.edge_grid → config.run → datastore).
        from .dataconvert import convert_data  # noqa: PLC0415

        source_collection_spec = dataspec[source_version][collection_name]
        data = read_data_collection(source_path, [source_collection_spec], file_codes)

        # Apply optional param overrides (inject missing Fortran log fields)
        if param_overrides:
            for key, value in param_overrides.items():
                for section in data["params"].values():
                    if isinstance(section, dict):
                        section[key] = value

        # Apply optional param aliases (rename Fortran parameter keys)
        if param_aliases:
            for py_name, fort_name in param_aliases.items():
                fort_lower = fort_name.lower()
                for section in data["params"].values():
                    if isinstance(section, dict) and fort_lower in section:
                        section[py_name] = section.pop(fort_lower)

        # Inject already-loaded micro_params so the converter has complete
        # parameters (particularly unit-bearing values).
        if self._micro_params is not None:
            if not isinstance(data["params"].get("micro_params"), dict):
                data["params"]["micro_params"] = {}
            data["params"]["micro_params"].update(self._micro_params.to_basedict())

        # ------------------------------------------------------------------
        # Convert to v2.0.0 (short-circuits if already at that version)
        # ------------------------------------------------------------------
        converted = convert_data(data, source_version, COMPATIBLE_DATASPEC_VERSION)

        # ------------------------------------------------------------------
        # Validate simulation count for per-simulation collections
        # ------------------------------------------------------------------
        if not target_spec.simulations_combined:
            for ds_name, ds_spec in target_spec.data.items():
                if ds_spec.data_location is None or ds_name not in converted:
                    continue
                n_source = len(converted[ds_name])
                n_target = self._macro_params.macro_simulations
                if n_source != n_target:
                    raise ValueError(
                        f"Source data has {n_source} simulation(s) but "
                        f"DataStore expects {n_target} "
                        f"(from macro_params.macro_simulations). "
                        "Simulation count must match."
                    )
                break  # one dataset is enough to confirm the count

        # ------------------------------------------------------------------
        # Write converted arrays into the HDF5, stamp provenance, and reload.
        #
        # Everything below mutates the file.  If any step fails, roll the
        # target collection back to its empty, freshly-initialized state so the
        # file never presents half-ingested data as "filled", then re-raise
        # loudly as ImportCollectionError.  Reads/conversion/validation above
        # do not touch the HDF5, so they need no rollback.
        # ------------------------------------------------------------------
        scale = "micro" if collection_name == "microscale_out" else "macro"
        try:
            for ds_name, ds_spec in target_spec.data.items():
                if ds_spec.data_location is None:
                    continue  # derived dataset — not stored on disk
                if ds_name not in converted:
                    continue  # optional dataset not produced by this converter
                if target_spec.simulations_combined:
                    self._write_array_to_dataset(ds_spec, converted[ds_name])
                else:
                    for sim_idx, arr in enumerate(converted[ds_name]):
                        self._write_array_to_dataset(ds_spec, arr, sim=sim_idx)

            # Record provenance if conversion happened
            if CONST.CONVERTED_FROM_ATTR in converted.get("params", {}):
                self._file.attrs[CONST.CONVERTED_FROM_ATTR] = (
                    converted["params"][CONST.CONVERTED_FROM_ATTR]
                )

            # Stamp pipeline provenance on the per-scale params group, plus
            # unconditional backend provenance (commit/dirty/compiler/type)
            # when an executable path is known.  Both stamps go onto the same
            # params group via the shared :meth:`stamp_provenance` helper.
            self.stamp_provenance(scale, "pipeline")
            if backend_executable is not None or replace_backend_attrs is not None:
                self.stamp_provenance(
                    scale,
                    "backend",
                    executable=backend_executable,
                    backend_override=backend_override,
                    replace_backend_attrs=replace_backend_attrs,
                )

            # Re-initialize in place (reloads all collections, params, etc.)
            self._file.flush()
            self._file.close()
            self.__init__(self._run_code, self._path, mode=self._mode)
        except Exception as exc:
            self._revert_collection_to_empty(collection_name)
            empty_state = (
                "MICRO_EMPTY" if scale == "micro" else "MACRO_EMPTY"
            )
            raise ImportCollectionError(
                f"Import of '{collection_name}' failed during the write phase "
                f"and was rolled back to {empty_state}; the HDF5 file holds no "
                f"partial data and the on-disk source is left uningested. "
                f"Original error: {exc}"
            ) from exc

    def _revert_collection_to_empty(self, collection_name):
        """Roll a partially-imported collection back to its empty state.

        Reuses the same machinery the ``init-*`` pathway uses, so the reverted
        file is equivalent to a freshly-initialized one — including the
        re-stamped ``init_*`` provenance:

        - ``microscale_out``: re-run :meth:`create` with ``force=True`` from
          the in-memory micro params.  No macroscale stage exists yet at
          micro-import time, so wiping and recreating the file is the correct
          reset.
        - ``macroscale_out``: re-run :meth:`initialize_macroscale` with
          ``force=True`` from the in-memory macro params, which wipes only the
          macroscale group and leaves the filled microscale data intact.

        The parameters come from the in-memory ``self._micro_params`` /
        ``self._macro_params`` (loaded from the HDF5 attributes at open time),
        which survive a failed write — no external file is needed.

        :param collection_name: ``"microscale_out"`` or ``"macroscale_out"``.
        :type collection_name: str
        """
        if collection_name == "microscale_out":
            micro = self._micro_params
            run_code, path, mode = self._run_code, self._path, self._mode
            if self._file:
                self._file.close()
            type(self).create(run_code, path, micro, force=True).close()
            self.__init__(run_code, path, mode=mode)
        else:  # macroscale_out
            self.initialize_macroscale(self._macro_params, force=True)

    # ------------------------------------------------------------------
    #  Derived collection helpers
    # ------------------------------------------------------------------

    def _can_generate_macroscale_in(self):
        """Check whether macroscale_in can be generated from available data.

        Requires microscale_out to be present, macro_params to be loaded,
        and the six datasets used by :func:`generate_macroscale_in` to exist
        with non-zero first dimension.

        :return: ``True`` if all preconditions are met.
        :rtype: bool
        """
        if "microscale_out" not in self._collections:
            return False
        if self._macro_params is None:
            return False

        required = [
            "pli_first_time",
            "tpa_leaving_time",
            "fiber_degraded",
            "sim_final_time",
            "tpa_unbound_by_pli",
            "tpa_unbound_kinetic",
        ]
        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
        for name in required:
            ds_spec = micro_spec.data.get(name)
            if ds_spec is None or ds_spec.data_location is None:
                return False
            if ds_spec.data_location not in self._file:
                return False
            if self._file[ds_spec.data_location].shape[0] == 0:
                return False
        return True

    def _build_macroscale_in_generator(self):
        """Build a closure that generates macroscale_in data on demand.

        The closure reads the six required microscale_out datasets from
        HDF5 (as numpy copies to avoid in-place mutation of HDF5 data),
        builds the params dict, and calls
        :func:`~lysis.dataio.dataconvert.generate_macroscale_in`.

        :return: A zero-argument callable returning a
            :class:`DataCollectionType` dict.
        :rtype: callable
        """
        # Capture references needed by the closure
        h5file = self._file
        micro_params = self._micro_params
        macro_params = self._macro_params
        micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]

        dataset_names = [
            "pli_first_time",
            "tpa_leaving_time",
            "fiber_degraded",
            "sim_final_time",
            "tpa_unbound_by_pli",
            "tpa_unbound_kinetic",
        ]

        def _generate():
            # Import inside closure to avoid circular imports
            from .dataconvert import generate_macroscale_in

            # Read datasets as numpy copies ([:] triggers full read)
            in_data = {}
            for name in dataset_names:
                ds_spec = micro_spec.data[name]
                in_data[name] = h5file[ds_spec.data_location][:]

            # Build params dict as expected by generate_macroscale_in
            in_data["params"] = {
                "micro_params": micro_params.to_basedict(),
                "macro_params": macro_params.to_basedict(),
            }

            return generate_macroscale_in(in_data)

        return _generate

    def close(self):
        """Close the underlying HDF5 file."""
        if self._file:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
