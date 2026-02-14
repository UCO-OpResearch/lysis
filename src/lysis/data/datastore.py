import os
import warnings

from enum import Flag, auto, unique
from typing import AnyStr

import h5py

from ..config.constants import CONST
from ..config.parameters import MicroParameters, MacroParameters
from .dataspec import DataCollectionSpec, DataSetSpec, dataspec
from .fileops import read_dataset

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


class DataStore:
    """Read-only interface to an HDF5 file following the v2.0.0 data specification.

    Provides dot-access to data collections and parameters::

        with DataStore("run01", "/path/to/data") as ds:
            ds.microscale_out.tpa_leaving_time[:10]  # lazy h5py.Dataset
            ds.macroscale_out[3].fiber_degrade_time   # per-simulation access
            ds.micro_params.fiber_radius              # loaded MicroParameters
            ds.macro_params.rows                      # loaded MacroParameters

    :param run_code: The run identifier (used to construct the HDF5 filename).
    :type run_code: str
    :param path: Directory containing the HDF5 file.
    :type path: str
    """

    def __init__(self, run_code, path):
        self._run_code = run_code
        self._path = path
        self._hdf5_path = os.path.join(path, f"{run_code}.h5")
        self._file = h5py.File(self._hdf5_path, "r")

        # Validate dataspec version
        found = self._file.attrs.get(CONST.DATASPEC_VERSION_ATTR)
        if found is None:
            self.close()
            raise ValueError(
                f"HDF5 file has no '{CONST.DATASPEC_VERSION_ATTR}' attribute: "
                f"{self._hdf5_path}"
            )
        if found != COMPATIBLE_DATASPEC_VERSION:
            self.close()
            raise ValueError(
                f"Dataspec version mismatch: file has '{found}', "
                f"expected '{COMPATIBLE_DATASPEC_VERSION}'"
            )
        self._dataspec_version = found

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

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in ("micro_params", "macro_params", "collections"):
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
