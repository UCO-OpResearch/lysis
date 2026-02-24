"""Data specification definitions for the lysis simulation system.

This module defines the schema for all data used in the lysis simulation,
including microscale and macroscale inputs and outputs. It provides:

- Type definitions for parameter and dataset structures
- Immutable specification classes (DataSetSpec, DataCollectionSpec)
- Complete data schemas for v1.99.0 (Fortran) and v2.0.0 (HDF5) formats
- Utility functions for shape parsing and data validation
- Tag-based version aliasing system

Data Organization
-----------------

The simulation has three main data stages:

1. **Microscale Output**: Results from fiber-level simulations
   - Stored per simulation run
   - Contains binding/unbinding statistics, lysis times, etc.

2. **Macroscale Input**: Processed microscale data for clot-level simulation
   - Binned degradation times and neighbor information
   - May be combined across multiple microscale runs

3. **Macroscale Output**: Results from clot-level simulations
   - Molecule locations, binding events, degradation state over time
   - Stored per simulation run with periodic snapshots

Specification Versions
----------------------

**v1.99.0 (Fortran format)**:
- File-based storage (text, binary, JSON)
- Separate files per dataset and simulation
- 1-based indexing for grid locations
- Legacy format for compatibility with Fortran code

**v2.0.0 (HDF5 format)**:
- Unified HDF5 container per run
- Attributes for parameters, datasets for data arrays
- 0-based indexing for grid locations
- Modern format with compression and metadata

Tag System
----------

Tags provide version aliasing for forward compatibility:
- ``"fortran"`` → ``"v1.99.0"``
- ``"hdf5"`` → ``"v2.0.0"``
- ``"current"`` → ``"hdf5"`` (can be changed as formats evolve)

Dynamic Shapes
--------------

Dataset shapes can reference parameter values using strings::

    shape=(-1, "macro_params.total_molecules")

At I/O time, parse_shape() resolves these to actual integers using
the parameter dictionary.

Usage Example
-------------

Accessing specifications::

    >>> from lysis.data.dataspec import dataspec
    >>> # Get v2.0.0 microscale output spec
    >>> micro_out_spec = dataspec["v2.0.0"]["microscale_out"]
    >>> # Access individual dataset spec
    >>> pli_spec = micro_out_spec.data["pli_first_time"]
    >>> pli_spec.dtype
    dtype('float64')
    >>> pli_spec.data_location
    'micro_data/pli_first_time'
    >>> # Use tag aliases
    >>> current_spec = dataspec["current"]["microscale_out"]

See Also
--------
- fileops.py : Functions that use these specs for actual I/O
- dataconvert.py : Functions that convert between spec versions
- datastore.py : High-level interface using these specs
"""

import copy
import dataclasses

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any, AnyStr, NewType, Union

import h5py
import numpy as np

from pint import Quantity

from ..config.constants import CONST, DataSetStorageType

__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# Type aliases for parameter dictionaries and datasets
BaseParamsType = NewType("BaseParamsType", dict[str, dict[str, int | float | str]])
"""Type for parameter dictionaries with base Python types only.

Structure: ``{"micro_params": {...}, "macro_params": {...}}``
where values are int, float, or str (no Pint Quantity objects).
Used for serialization to JSON and other formats that don't support custom types.
"""

UnitParamsType = NewType(
    "UnitParamsType", dict[str, dict[str, int | float | str | Quantity]]
)
"""Type for parameter dictionaries that may include Pint Quantity objects.

Structure: ``{"micro_params": {...}, "macro_params": {...}}``
where values can be int, float, str, or Pint Quantity with units.
Used internally when working with parameters that have physical units.
"""

DataSetType = Union[np.ndarray | list[np.ndarray] | BaseParamsType]
"""Type for individual datasets.

Can be:
- Single numpy array (for combined simulations)
- List of numpy arrays (for separate simulations)
- Parameter dictionary (for parameter datasets)
"""

DataCollectionType = NewType("DataCollectionType", dict[str, DataSetType])
"""Type for a collection of datasets.

Dictionary mapping dataset names to their data (arrays or parameter dicts).
Represents all data in a collection (e.g., all microscale output datasets).
"""


@dataclass(frozen=True, kw_only=True)
class DataSetSpec:
    """Specification for a single dataset's storage and structure.

    Defines how a dataset is stored (format, location) and what it contains
    (data type, shape). Frozen to ensure specs don't change at runtime.

    :ivar dataset_storage_type: Storage backend for this dataset
    :vartype dataset_storage_type: DataSetStorageType
    :ivar dtype: NumPy dtype or special type (e.g., h5py.string_dtype(), Quantity)
    :vartype dtype: np.dtype
    :ivar data_location: Path template for file-based storage or HDF5 group/dataset path.
        May contain format placeholders like {sim:02} or {file_code}.
        None for derived datasets not stored directly.
    :vartype data_location: str | None
    :ivar shape: Expected array shape. Tuple of integers or strings.
        -1 indicates variable dimension.
        Strings reference parameter values (e.g., "macro_params.total_molecules").
        Default (-1,) means 1D array with variable length.
    :vartype shape: tuple[int | str, ...]
    :ivar delimiter: Delimiter for text files (e.g., "," for CSV, "\\u0000" for null-terminated).
        None for binary/HDF5 formats.
    :vartype delimiter: str | None

    Examples
    --------
    Fixed-size array stored in HDF5::

        DataSetSpec(
            dataset_storage_type=DataSetStorageType.HDF5_DATASET,
            data_location="micro_data/pli_first_time",
            dtype=np.float64,
            shape=()  # Scalar
        )

    Variable-size array with dynamic dimension from parameters::

        DataSetSpec(
            dataset_storage_type=DataSetStorageType.FILE_BINARY,
            data_location="{sim:02}/m_loc{file_code}_{sim:02}.dat",
            dtype=np.int32,
            shape=(-1, "macro_params.total_molecules")  # Rows variable, cols from params
        )

    Structured array (event log)::

        DataSetSpec(
            dataset_storage_type=DataSetStorageType.FILE_TEXT,
            data_location="{sim:02}/f_deg_list{file_code}_{sim:02}.dat",
            dtype=np.dtype([
                ("Simulation Time Elapsed", np.float64),
                ("Grid Location Index", np.int32),
                ("Fiber New Degrade Time", np.float64)
            ]),
            delimiter=","
        )
    """

    dataset_storage_type: DataSetStorageType
    dtype: np.dtype
    data_location: str | None = None
    shape: tuple[int, ...] = (-1,)
    delimiter: str | None = None
    # Hidden fields — not in constructor, auto-populated by DataSpec
    version: str = field(init=False, default="")
    collection: str = field(init=False, default="")
    name: str = field(init=False, default="")

    def __copy__(self):
        # dataclasses.replace creates a new instance with the same init field values.
        # Hidden fields (init=False) must be copied manually since they're not in __init__.
        new = dataclasses.replace(self)
        object.__setattr__(new, "version", self.version)
        object.__setattr__(new, "collection", self.collection)
        object.__setattr__(new, "name", self.name)
        return new

    def __deepcopy__(self, memo):
        new = dataclasses.replace(
            self,
            dtype=copy.deepcopy(self.dtype, memo),
            shape=copy.deepcopy(self.shape, memo),
        )
        object.__setattr__(new, "version", self.version)
        object.__setattr__(new, "collection", self.collection)
        object.__setattr__(new, "name", self.name)
        memo[id(self)] = new
        return new


@dataclass(frozen=True, kw_only=True)
class DataCollectionSpec:
    """Specification for a collection of related datasets.

    Groups multiple datasets that are logically related (e.g., all microscale
    outputs, all macroscale inputs). Defines whether simulations are stored
    together or separately. Frozen to ensure specs don't change at runtime.

    :ivar simulations_combined: Whether multiple simulations are in one file/group.
        - True: All simulation data combined in single location (e.g., v1.99.0 microscale output)
        - False: Each simulation in separate location (e.g., v1.99.0 macroscale output with {sim:02})
    :vartype simulations_combined: bool
    :ivar params: Specification for parameter storage in this collection.
        None if collection has no associated parameters (e.g., macroscale input
        which uses parameters from microscale output).
    :vartype params: DataSetSpec | None
    :ivar data: Dictionary mapping dataset names to their specifications.
        Keys are dataset names (e.g., "pli_first_time", "snapshot_time").
        Values are DataSetSpec objects defining storage and structure.
    :vartype data: dict[str, DataSetSpec]

    Examples
    --------
    Microscale output (combined simulations, has parameters)::

        DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=DataSetStorageType.FILE_JSON,
                dtype=Quantity
            ),
            data={
                "pli_first_time": DataSetSpec(...),
                "fiber_degraded": DataSetSpec(...),
                # ... more datasets
            }
        )

    Macroscale output (separate simulations, has parameters)::

        DataCollectionSpec(
            simulations_combined=False,  # Separate files per simulation
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=DataSetStorageType.FILE_JSON,
                dtype=Quantity
            ),
            data={
                "snapshot_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/snapshot_time",  # Note {sim:02}
                    ...
                ),
                # ... more datasets
            }
        )

    Macroscale input (combined, no parameters of its own)::

        DataCollectionSpec(
            simulations_combined=True,
            params=None,  # Uses parameters from microscale output
            data={
                "bin_edge_proportions": DataSetSpec(...),
                # ... more datasets
            }
        )
    """

    simulations_combined: bool
    params: DataSetSpec
    data: dict[str, DataSetSpec]
    # Hidden fields — not in constructor, auto-populated by DataSpec
    version: str = field(init=False, default="")
    collection: str = field(init=False, default="")

    def replace(self, **changes):
        """Create a modified copy of this collection spec.

        A recursive version of :func:`dataclasses.replace` that accepts changes
        both to the fields of this :class:`DataCollectionSpec` *and* to the
        :class:`DataSetSpec` entries stored inside ``self.data``. The original
        object is never mutated.

        Each keyword argument key must be one of:

        - A named *init* field of :class:`DataCollectionSpec`
          (``simulations_combined``, ``params``, or ``data``).
        - A key present in ``self.data`` (i.e. a dataset name such as
          ``"pli_first_time"``).

        Each keyword argument value may be one of:

        - A replacement object used directly (a :class:`DataSetSpec` or a
          :class:`bool` for boolean fields such as ``simulations_combined``).
        - A :class:`dict` of field changes forwarded to
          :func:`dataclasses.replace` on the *current* child
          :class:`DataSetSpec` — unspecified fields keep their original values.

        :param changes: Keyword arguments mapping field names or dataset names
            to their new values or dicts of sub-field changes.
        :raises KeyError: If a key is not an init field of this class and not a
            key in ``self.data``.
        :raises ValueError: If a value is not a :class:`DataSetSpec`, :class:`bool`,
            or :class:`dict`.
        :return: A new, frozen :class:`DataCollectionSpec` with the specified
            changes applied. All unchanged ``data`` entries are deep-copied.
        :rtype: DataCollectionSpec

        Examples
        --------
        Toggle the ``simulations_combined`` flag::

            new_coll = coll.replace(simulations_combined=False)

        Replace the ``params`` spec with a new :class:`DataSetSpec`::

            new_coll = coll.replace(
                params=DataSetSpec(
                    dataset_storage_type=DataSetStorageType.FILE_JSON,
                    dtype=float,
                )
            )

        Update only selected fields of the existing ``params`` spec::

            new_coll = coll.replace(params={"dtype": np.float64})

        Replace a named dataset entry with a new :class:`DataSetSpec`::

            new_coll = coll.replace(pli_first_time=DataSetSpec(...))

        Update only selected fields of a named dataset entry::

            new_coll = coll.replace(pli_first_time={"dtype": np.float64})

        Apply multiple changes in one call::

            new_coll = coll.replace(
                simulations_combined=False,
                params={"dtype": np.float64},
                pli_first_time={"data_location": "micro_data/pli_first_time"},
            )

        Replace the entire ``data`` dict (all entries at once)::

            new_coll = coll.replace(
                data={"pli_first_time": DataSetSpec(...), "sim_final_time": DataSetSpec(...)}
            )

        """
        own_changes = {}
        if "data" in changes.keys():
            # Caller supplied a full replacement dict: deep-copy it so the new
            # spec is fully independent of the caller's object, then remove the
            # key from changes so the loop below does not encounter "data"
            # (whose generic type dict[str, DataSetSpec] is not compatible with
            # isinstance()).
            new_data = copy.deepcopy(changes.pop("data"))
        else:
            # No full replacement supplied: start from a deep copy of the current
            # data dict so that individual entry changes below are isolated from
            # the original and unchanged entries are not shared with the new spec.
            new_data = copy.deepcopy(self.data)

        # Build a {field_name: field_type} lookup for every constructor-visible
        # field.  Using the declared type lets us validate values via
        # isinstance(v, field_types[k]) without a per-field if-statement.
        #
        # Note: the "data" field has type dict[str, DataSetSpec], a generic alias
        # that isinstance() cannot accept.  That field is always popped from
        # `changes` before this loop, so field_types["data"] is never passed to
        # isinstance() here.
        field_types = {f.name: f.type for f in dataclasses.fields(self) if f.init}

        for k, v in changes.items():
            if k in field_types:
                # --- Change targets a field on this DataCollectionSpec ---
                if isinstance(v, dict):
                    # Dict form: apply the sub-field changes to the existing child
                    # DataSetSpec via dataclasses.replace, preserving all other fields.
                    own_changes[k] = dataclasses.replace(getattr(self, k), **v)
                elif isinstance(v, field_types[k]):
                    # Direct replacement: value must match the field's declared type.
                    own_changes[k] = v
                else:
                    raise ValueError(
                        f"Expected a {field_types[k].__name__} or dict for field "
                        f"'{k}', got {type(v).__name__}: {v!r}"
                    )
            elif k in self.data.keys():
                # --- Change targets a named dataset inside self.data ---
                if isinstance(v, DataSetSpec):
                    # Use the provided value directly, replacing the existing entry.
                    new_data[k] = v
                elif isinstance(v, dict):
                    # Dict form: apply the sub-field changes to the existing entry,
                    # based on the *original* spec (self.data[k]), not the copy.
                    new_data[k] = dataclasses.replace(self.data[k], **v)
                else:
                    raise ValueError(
                        f"Changes should be fields of this class, this class's data, "
                        f"or a dict of changes for a DataSetSpec. {k, v}"
                    )
            else:
                raise KeyError(k)

        # Construct the new frozen DataCollectionSpec with the updated fields and
        # data dict; dataclasses.replace handles the frozen constraint via object.__setattr__.
        return dataclasses.replace(self, data=new_data, **own_changes)

    def __copy__(self):
        # dataclasses.replace creates a new instance with the same init field values.
        # Hidden fields (init=False) must be copied manually since they're not in __init__.
        new = dataclasses.replace(self)
        object.__setattr__(new, "version", self.version)
        object.__setattr__(new, "collection", self.collection)
        return new

    def __deepcopy__(self, memo):
        new = dataclasses.replace(
            self,
            params=copy.deepcopy(self.params, memo),
            data={k: copy.deepcopy(v, memo) for k, v in self.data.items()},
        )
        object.__setattr__(new, "version", self.version)
        object.__setattr__(new, "collection", self.collection)
        memo[id(self)] = new
        return new


class DataSpec:
    """Versioned specification containing a set of data collections.

    Wraps a dictionary of :class:`DataCollectionSpec` objects and carries the
    version string. Auto-populates hidden fields (``.version``, ``.collection``,
    ``.name``) on child specs so they know their own identity.

    Supports dict-like access for backward compatibility::

        spec = DataSpec("v2.0.0", {"microscale_out": ..., ...})
        spec["microscale_out"]         # DataCollectionSpec
        spec.version                   # "v2.0.0"
        for name, coll in spec.items(): ...

    :param version: The version string (e.g., ``"v2.0.0"``).
    :type version: str
    :param collections: Mapping of collection names to their specifications.
    :type collections: dict[str, DataCollectionSpec]
    """

    def __init__(self, version: str, collections: dict[str, DataCollectionSpec]):
        self._version = version
        self._collections: dict[str, DataCollectionSpec] = {}
        for coll_name, coll_spec in collections.items():
            # Populate hidden fields on the frozen DataCollectionSpec
            object.__setattr__(coll_spec, "version", version)
            object.__setattr__(coll_spec, "collection", coll_name)
            # Populate hidden fields on the params DataSetSpec
            if coll_spec.params is not None:
                object.__setattr__(coll_spec.params, "version", version)
                object.__setattr__(coll_spec.params, "collection", coll_name)
                object.__setattr__(coll_spec.params, "name", "params")
            # Populate hidden fields on each dataset DataSetSpec
            for ds_name, ds_spec in coll_spec.data.items():
                object.__setattr__(ds_spec, "version", version)
                object.__setattr__(ds_spec, "collection", coll_name)
                object.__setattr__(ds_spec, "name", ds_name)
            self._collections[coll_name] = coll_spec

    @property
    def version(self) -> str:
        """The version string for this specification."""
        return self._version

    def __getitem__(self, key: str) -> DataCollectionSpec:
        return self._collections[key]

    def __contains__(self, key: str) -> bool:
        return key in self._collections

    def __iter__(self):
        return iter(self._collections)

    def __len__(self) -> int:
        return len(self._collections)

    def items(self):
        return self._collections.items()

    def keys(self):
        return self._collections.keys()

    def values(self):
        return self._collections.values()

    def __str__(self) -> str:
        colls = ", ".join(self._collections.keys())
        return f"<DataSpec '{self._version}', collections=[{colls}]>"

    def __repr__(self) -> str:
        colls = ", ".join(c.__repr__() for c in self._collections.values())
        return f"<DataSpec '{self._version}', collections=[{colls}]>"


# Master data specification dictionary
# Structure: dataspec[version][collection_name] -> DataCollectionSpec
# Versions: "v1.99.0" (Fortran), "v2.0.0" (HDF5), plus tag aliases
# Collections: "microscale_out", "macroscale_in", "macroscale_out"
_dataspec_raw: dict[str, dict[str, DataCollectionSpec]] = {
    # ============================================================================
    # v1.95.0: Identical to v1.99.0, except pre Pint Quantity implementation
    # NOTE: This spec will be added below as a copy of v1.99.0
    # ============================================================================
    "v1.95.0": {},
    # ============================================================================
    # v1.99.0: Fortran-compatible file-based format
    # ============================================================================
    "v1.99.0": {
        # ------------------------------------------------------------------------
        # Microscale output: Results from fiber-level simulations
        # All simulations combined in single files
        # ------------------------------------------------------------------------
        "microscale_out": DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                # Log file with simulation status messages
                "micro_log": DataSetSpec(
                    data_location="micro{file_code}.txt",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",  # Null-terminated strings
                ),
                # Time when first plasmin (PLi) molecule was generated in each simulation
                "firstPLi": DataSetSpec(
                    data_location="firstPLi{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                # Number of tPA molecules still in fiber at simulation end
                "lasttPA": DataSetSpec(
                    data_location="lasttPA{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                # Whether lysis completed (>= snap_proportion of doublets degraded)
                "lyscomplete": DataSetSpec(
                    data_location="lyscomplete{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.uint32,
                ),
                # Time elapsed in each simulation (lysis time if completed)
                "lysis": DataSetSpec(
                    data_location="lysis{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                # Number of plasmin molecules generated by simulation end
                "PLi": DataSetSpec(
                    data_location="PLi{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                # Time when tPA leaves system (or infinity if still bound)
                "tPA_time": DataSetSpec(
                    data_location="tPA_time{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                # Whether tPA was micro-unbound (forced unbinding by plasmin degradation)
                "tPAPLiunbd": DataSetSpec(
                    data_location="tPAPLiunbd{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                # Whether tPA kinetically unbound (not forced, free to rebind)
                "tPAunbind": DataSetSpec(
                    data_location="tPAunbind{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
            },
        ),
        # ------------------------------------------------------------------------
        # Macroscale input: Processed microscale data for clot-level simulation
        # Binned degradation times and neighbor information
        # All simulations combined; no separate parameters (uses microscale params)
        # ------------------------------------------------------------------------
        "macroscale_in": DataCollectionSpec(
            simulations_combined=True,
            params=None,  # Uses parameters from microscale output
            data={
                # Bin boundaries of CDF for tPA leaving time distribution (100 bins, 101 edges)
                "tPAleave": DataSetSpec(
                    data_location="tPAleave{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                # Times at bin boundaries: i% of simulations have tPA leave by tsectPA[i]
                "tsectPA": DataSetSpec(
                    data_location="tsectPA{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                # Fiber degrade times, binned by tPA leaving time, sorted within each bin
                "lysismat": DataSetSpec(
                    data_location="lysismat{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(-1, 100),
                ),
                # Count of simulations per bin where full fiber lysis occurred
                "lenlysisvect": DataSetSpec(
                    data_location="lenlysisvect{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(100,),
                ),
                # 1-indexed Fortran locations of neighboring fibers for each edge in grid
                "neighbors": DataSetSpec(
                    data_location="neighbors{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.int32,
                    shape=(-1,),
                ),
            },
        ),
        # ------------------------------------------------------------------------
        # Macroscale output: Results from clot-level simulations
        # Each simulation stored in separate files (note {sim:02} in paths)
        # Contains snapshots of molecule locations, binding events, degradation state
        # ------------------------------------------------------------------------
        "macroscale_out": DataCollectionSpec(
            simulations_combined=False,  # Separate files per simulation
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                # Log file with macroscale simulation status messages
                "macro_log": DataSetSpec(
                    data_location="{sim:02}/macro{file_code}_{sim:02}.txt",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                # Number of snapshots recorded during this simulation
                "Nsave": DataSetSpec(
                    data_location="{sim:02}/Nsave{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(),
                ),
                # Simulation time at each snapshot
                "tsave": DataSetSpec(
                    data_location="{sim:02}/tsave{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                # Event log: time, location, new degrade time for each fiber degradation update
                "f_deg_list": DataSetSpec(
                    data_location="{sim:02}/f_deg_list{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("Grid Location Index", np.int32),
                            ("Fiber New Degrade Time", np.float64),
                        ]
                    ),
                    delimiter=",",
                ),
                # Event log: time, molecule ID, new MolStatus, location for each binding/unbinding event
                "m_bind_t": DataSetSpec(
                    data_location="{sim:02}/m_bind_t{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("tPA Molecule Index", np.int32),
                            ("Molecule New Status", np.int32),
                            ("Grid Location Index", np.int32),
                        ]
                    ),
                    delimiter=",",
                ),
                # Fortran location (1-indexed) of each tPA molecule at each snapshot
                "m_loc": DataSetSpec(
                    data_location="{sim:02}/m_loc{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(-1, "macro_params.total_molecules"),
                ),
                # Bound/unbound status (1/0) of each tPA molecule at each snapshot
                "m_bound": DataSetSpec(
                    data_location="{sim:02}/m_bound{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(-1, "macro_params.total_molecules"),
                ),
                # Time when each tPA molecule first reached the back row (first passage times)
                "mfpt": DataSetSpec(
                    data_location="{sim:02}/mfpt{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
            },
        ),
    },
    # ============================================================================
    # v2.0.0: Modern HDF5-based unified storage format
    # ============================================================================
    "v2.0.0": {
        # ------------------------------------------------------------------------
        # Microscale output: Results from fiber-level simulations in HDF5
        # All simulations combined in single HDF5 file
        # Parameters stored as HDF5 attributes on micro_data group
        # ------------------------------------------------------------------------
        "microscale_out": DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="micro_data",  # HDF5 group path
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                # Log file with microscale simulation status messages
                "micro_log": DataSetSpec(
                    data_location="log_files/micro_log",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=h5py.string_dtype(),
                ),
                # Time when first plasmin (PLi) molecule was generated in each simulation
                "pli_first_time": DataSetSpec(
                    data_location="micro_data/pli_first_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                # Number of tPA molecules still in fiber at simulation end
                "tpa_final_num": DataSetSpec(
                    data_location="micro_data/tpa_final_num",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.uint8,
                ),
                # Whether lysis completed (>= snap_proportion of doublets degraded)
                "fiber_degraded": DataSetSpec(
                    data_location="micro_data/fiber_degraded",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                # Time elapsed in each simulation (lysis time if completed)
                "sim_final_time": DataSetSpec(
                    data_location="micro_data/sim_final_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                # Number of plasmin molecules generated by simulation end
                "pli_generated_num": DataSetSpec(
                    data_location="micro_data/pli_generated_num",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.uint16,
                ),
                # Time when tPA leaves system (or infinity if still bound)
                "tpa_leaving_time": DataSetSpec(
                    data_location="micro_data/tpa_leaving_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                # Whether tPA was micro-unbound (forced unbinding by plasmin degradation)
                "tpa_unbound_by_pli": DataSetSpec(
                    data_location="micro_data/tpa_unbound_by_pli",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                # Whether tPA kinetically unbound (not forced, free to rebind)
                "tpa_unbound_kinetic": DataSetSpec(
                    data_location="micro_data/tpa_unbound_kinetic",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
            },
        ),
        # ------------------------------------------------------------------------
        # Macroscale input: Processed microscale data for clot-level simulation
        # These are computed/derived datasets, not stored on disk (data_location=None)
        # Generated by dataconvert.generate_macroscale_in() from microscale output
        # ------------------------------------------------------------------------
        "macroscale_in": DataCollectionSpec(
            simulations_combined=True,
            params=None,  # Uses parameters from microscale output
            data={
                # Bin boundaries of CDF for tPA leaving time distribution (100 bins, 101 edges)
                "bin_edge_proportions": DataSetSpec(
                    data_location=None,  # Derived dataset, not stored
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(101,),
                ),
                # Times at bin boundaries: i% of simulations have tPA leave by this time
                "bin_edge_tpa_leaving_time": DataSetSpec(
                    data_location=None,  # Derived dataset, not stored
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(101,),
                ),
                # Fiber degrade times, binned by tPA leaving time, sorted within each bin
                "binned_fiber_degrade_time": DataSetSpec(
                    data_location=None,  # Derived dataset, not stored
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(-1, 100),
                ),
                # Count of simulations per bin where full fiber lysis occurred
                "binned_fiber_degraded": DataSetSpec(
                    data_location=None,  # Derived dataset, not stored
                    dataset_storage_type=None,
                    dtype=np.uint16,
                    shape=(100,),
                ),
                # 0-indexed locations of neighboring edges for each edge in grid (8 neighbors each)
                "edge_grid_neighbors": DataSetSpec(
                    data_location=None,  # Derived dataset, not stored
                    dataset_storage_type=None,
                    dtype=np.uint32,
                    shape=(-1, 8),
                ),
            },
        ),
        # ------------------------------------------------------------------------
        # Macroscale output: Results from clot-level simulations in HDF5
        # Each simulation stored in separate HDF5 groups (sim_{sim:02})
        # Parameters stored as attributes on macro_data group
        # Uses 2D grid indexing (row, rank) instead of 1D Fortran indices
        # ------------------------------------------------------------------------
        "macroscale_out": DataCollectionSpec(
            simulations_combined=False,  # Separate groups per simulation
            params=DataSetSpec(
                data_location="macro_data",  # HDF5 group path
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                # Log file with macroscale simulation status messages
                "macro_log": DataSetSpec(
                    data_location="log_files/macro_log__sim_{sim:02}",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=h5py.string_dtype(),
                ),
                # Simulation time at each snapshot
                "snapshot_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/snapshot_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                # Event log: time, location (row, rank), new degrade time for each fiber update
                "fiber_degrade_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/fiber_degrade_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("Grid Location Row", np.uint32),
                            ("Grid Location Rank", np.uint32),
                            ("Fiber New Degrade Time", np.float64),
                        ]
                    ),
                ),
                # Event log: time, molecule ID, new MolStatus, location (row, rank) for each binding/unbinding
                "tpa_bind_events": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_bind_events",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("tPA Molecule Index", np.uint64),
                            (
                                "Molecule New Status",
                                h5py.enum_dtype(
                                    {i.name: i.value for i in CONST.MOL_STATUS},
                                    basetype="u1",
                                ),
                            ),
                            ("Grid Location Row", np.uint32),
                            ("Grid Location Rank", np.uint32),
                        ]
                    ),
                ),
                # Grid location (row, rank) of each tPA molecule at each snapshot
                "tpa_location_snapshot": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_location_snapshot",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.uint32,
                    shape=(-1, 2, -1),
                ),
                # Time when each tPA molecule first reached the back row (first passage times)
                "tpa_transit_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_transit_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
            },
        ),
    },
}


# ============================================================================
# Adding Derived DataSpecs
# ============================================================================
# v1.95.00 <-- v1.99.0
def _create_v1_95():
    v1_95 = copy.deepcopy(_dataspec_raw["v1.99.0"])
    v1_95["microscale_out"] = v1_95["microscale_out"].replace(
        params=DataSetSpec(
            data_location="micro{file_code}.txt",
            dtype=np.float64,
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_PARSED,
        )
    )
    v1_95["macroscale_out"] = v1_95["macroscale_out"].replace(
        params={"dtype": np.float64}
    )
    return v1_95


_dataspec_raw["v1.95.0"] = _create_v1_95()


# Wrap each version's raw dict in a DataSpec object
dataspec: dict[str, DataSpec] = {}
for _version, _collections in _dataspec_raw.items():
    dataspec[_version] = DataSpec(_version, _collections)
del _version, _collections

# ============================================================================
# Tag System: Version Aliases for Forward Compatibility
# ============================================================================
# Tags allow code to reference versions without hard-coding version numbers.
# For example, code can use dataspec["current"] and when the format evolves,
# only this dictionary needs to be updated.
tags = {
    "fortran": "v1.99.0",  # Alias for Fortran-compatible file format
    "hdf5": "v2.0.0",  # Alias for HDF5 unified storage format
    "current": "hdf5",  # Alias for current recommended format (can be changed)
}

# Add tag aliases to the dataspec dictionary by copying the referenced DataSpec objects
for _k, _v in tags.items():
    dataspec[_k] = dataspec[_v]
del _k, _v

# Spec versions that use Fortran file-based storage.
# Derived from the tag system so this set stays in sync if tags are updated.
# Used by fileops.read_data_collection to trigger paramcheck validation at read time.
fortran_versions: frozenset[str] = frozenset({tags["fortran"], "v1.95.0"})


def parse_shape(
    shape: tuple[int | str, ...], params: BaseParamsType = None
) -> tuple[int, ...]:
    """Resolve dynamic shape dimensions using parameter values.

    Converts shape specifications that may contain parameter references (strings)
    into concrete integer tuples suitable for NumPy array operations. This allows
    dataset shapes to depend on simulation parameters like molecule counts or
    grid dimensions.

    Shape elements can be:
    - **Integers**: Passed through unchanged (e.g., 100, -1)
    - **Strings**: Parameter references in "dict_key.param_name" format
      (e.g., "macro_params.total_molecules")

    :param shape: Shape tuple from DataSetSpec, may contain integers and/or strings.
        Integers (including -1 for variable dimensions) are kept as-is.
        Strings are parameter references resolved via params dictionary.
    :type shape: tuple[int | str, ...]
    :param params: Parameter dictionary with structure::

            {
                "micro_params": {"param1": value1, ...},
                "macro_params": {"param2": value2, ...}
            }

        Required if shape contains any string references.
    :type params: BaseParamsType, optional
    :raises RuntimeError: If shape contains a type other than int or str
    :return: Fully resolved shape tuple with all integer dimensions
    :rtype: tuple[int, ...]

    Examples
    --------
    Fixed shape (no parameters needed)::

        >>> parse_shape((100, 3))
        (100, 3)

    Variable first dimension::

        >>> parse_shape((-1, 100))
        (-1, 100)

    Shape with parameter reference::

        >>> params = {"macro_params": {"total_molecules": 43074}}
        >>> parse_shape((-1, "macro_params.total_molecules"), params)
        (-1, 43074)

    Multiple parameter references::

        >>> params = {
        ...     "macro_params": {"total_molecules": 43074, "number_of_saves": 121}
        ... }
        >>> parse_shape(("macro_params.number_of_saves", "macro_params.total_molecules"), params)
        (121, 43074)

    See Also
    --------
    check_dataset_spec : Validates data against spec using parsed shapes
    DataSetSpec : Contains shape specifications
    """
    parsed_shape = []
    for i in shape:
        if isinstance(i, int):
            # Integer dimension (fixed size or -1 for variable)
            parsed_shape.append(i)
        elif isinstance(i, str):
            # String dimension - resolve from parameter dictionary
            # Format: "dict_name.param_name" e.g., "macro_params.total_molecules"
            parts = i.split(".")
            parsed_shape.append(params[parts[0]][parts[1]])
        else:
            raise RuntimeError(f"Incorrect shape format: {i}. Expected int or str.")
    return tuple(parsed_shape)


def check_dataset_spec(
    data: np.ndarray, spec: DataSetSpec, params: BaseParamsType = None
) -> bool:
    """Validate that a data array conforms to its specification.

    Checks two aspects of conformance:
    1. **Data type compatibility**: Can data.dtype be cast to spec.dtype?
    2. **Shape matching**: Does data.shape match spec.shape (after resolving parameters)?

    Shape validation allows -1 in spec.shape to indicate variable dimensions
    that will match any size in the actual data.

    :param data: The data array to validate
    :type data: np.ndarray
    :param spec: The specification defining expected dtype and shape
    :type spec: DataSetSpec
    :param params: Parameter dictionary for resolving dynamic shape dimensions.
        Required if spec.shape contains string references.
    :type params: BaseParamsType, optional
    :return: True if data conforms to spec (dtype and shape match), False otherwise
    :rtype: bool

    Examples
    --------
    Fixed shape validation::

        >>> import numpy as np
        >>> spec = DataSetSpec(
        ...     dataset_storage_type=DataSetStorageType.HDF5_DATASET,
        ...     dtype=np.float64,
        ...     shape=(100, 3)
        ... )
        >>> data = np.zeros((100, 3))
        >>> check_dataset_spec(data, spec)
        True
        >>> bad_data = np.zeros((100, 4))  # Wrong shape
        >>> check_dataset_spec(bad_data, spec)
        False

    Variable dimension (any size accepted)::

        >>> spec = DataSetSpec(
        ...     dataset_storage_type=DataSetStorageType.HDF5_DATASET,
        ...     dtype=np.int32,
        ...     shape=(-1, 100)  # First dimension variable
        ... )
        >>> check_dataset_spec(np.zeros((50, 100), dtype=np.int32), spec)
        True
        >>> check_dataset_spec(np.zeros((1000, 100), dtype=np.int32), spec)
        True
        >>> check_dataset_spec(np.zeros((100, 50), dtype=np.int32), spec)
        False  # Second dimension doesn't match

    Dynamic shape from parameters::

        >>> params = {"macro_params": {"total_molecules": 43074}}
        >>> spec = DataSetSpec(
        ...     dataset_storage_type=DataSetStorageType.FILE_BINARY,
        ...     dtype=np.int32,
        ...     shape=(-1, "macro_params.total_molecules")
        ... )
        >>> data = np.zeros((121, 43074), dtype=np.int32)
        >>> check_dataset_spec(data, spec, params)
        True

    Data type compatibility::

        >>> spec = DataSetSpec(
        ...     dataset_storage_type=DataSetStorageType.HDF5_DATASET,
        ...     dtype=np.float64,
        ...     shape=(100,)
        ... )
        >>> check_dataset_spec(np.zeros(100, dtype=np.float32), spec)
        True  # float32 can be cast to float64
        >>> check_dataset_spec(np.zeros(100, dtype=np.int32), spec)
        True  # int32 can be cast to float64
        >>> spec_int = DataSetSpec(
        ...     dataset_storage_type=DataSetStorageType.HDF5_DATASET,
        ...     dtype=np.int32,
        ...     shape=(100,)
        ... )
        >>> check_dataset_spec(np.zeros(100, dtype=np.float64), spec_int)
        False  # float64 cannot safely cast to int32

    See Also
    --------
    parse_shape : Resolves dynamic dimensions in shape specifications
    DataSetSpec : Specification class containing dtype and shape requirements
    """
    # Check data type compatibility
    if not np.can_cast(data.dtype, spec.dtype):
        return False

    # Check shape compatibility (resolve any parameter references first)
    for idx, expected_dim in enumerate(parse_shape(spec.shape, params=params)):
        if expected_dim < 0:
            # -1 means variable dimension, accept any size
            continue
        elif expected_dim != data.shape[idx]:
            # Fixed dimension must match exactly
            return False

    return True
