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

from enum import Flag, auto, unique
from typing import AnyStr

import h5py

from ..config.constants import CONST
from ..config.parameters import MicroParameters, MacroParameters
from .dataspec import DataCollectionSpec, DataSetSpec, dataspec, parse_shape
from .fileops import read_dataset, write_dataset, _validate_hdf5_version

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

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


@unique
class DataStatus(Flag):
    NONE = 0
    INITIALIZED = auto()
    LOADED = auto()
    SAVED = auto()
    FILLED = auto()


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
    """Read-only view of one simulation's datasets within a per-sim collection.

    Provides dot-access to HDF5 datasets for a single simulation index.
    Dataset names are resolved via the spec's ``data_location`` field,
    formatted with the simulation index.

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
    """Read-only interface to one data collection within the HDF5 file.

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
    the supplied *generator* callable, and cached for subsequent accesses.

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
        """Call the generator if data has not yet been generated."""
        if self._data is None:
            self._data = self._generator()

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

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in ("mode", "micro_params", "macro_params", "collections", "status"):
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
    def create(cls, run_code, path, micro_params):
        """Create a new DataStore with microscale parameters and empty datasets.

        Creates a new HDF5 file, writes the dataspec version attribute,
        stores ``micro_params`` as HDF5 attributes on the ``micro_data``
        group, and creates zero-length datasets for all microscale_out
        datasets defined in the v2.0.0 specification.

        :param run_code: The Run identifier (used to construct the HDF5
            filename as ``{run_code}.h5``).
        :type run_code: str
        :param path: Directory where the HDF5 file will be created.
        :type path: str
        :param micro_params: The microscale parameters to store.
        :type micro_params: MicroParameters
        :return: A new DataStore opened in ``"a"`` (read/write) mode with
            microscale_out collection.
        :rtype: DataStore
        :raises FileExistsError: If the HDF5 file already exists.
        :raises TypeError: If ``micro_params`` is not a
            :class:`~lysis.config.parameters.MicroParameters` instance.
        """
        if not isinstance(micro_params, MicroParameters):
            raise TypeError(
                f"Expected MicroParameters, got {type(micro_params).__name__}"
            )

        hdf5_path = os.path.join(path, f"{run_code}.h5")
        if os.path.exists(hdf5_path):
            raise FileExistsError(f"HDF5 file already exists: {hdf5_path}")

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

        return cls(run_code, path, mode="a")

    def initialize_macroscale(self, macro_params):
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

        Modifies the DataStore **in place** and returns ``None``.
        Requires the DataStore to be opened in a writable mode
        (e.g., ``"a"``)::

            with DataStore(run_code, path, mode="a") as ds:
                ds.initialize_macroscale(run.macro_params)

        :param macro_params: The macroscale parameters to store.  All
            fields except ``forced_unbind`` are used as provided.
        :type macro_params: MacroParameters
        :raises IOError: If the DataStore is opened in read-only mode.
        :raises TypeError: If ``macro_params`` is not a
            :class:`~lysis.config.parameters.MacroParameters` instance.
        :raises ValueError: If microscale_out is not present, if
            macroscale_out is already present, or if the microscale
            datasets are empty (i.e. no Simulations have been written).
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
            raise ValueError(
                "macroscale_out is already present in this DataStore."
            )

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
