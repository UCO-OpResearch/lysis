"""Constants, enumerations, and configuration values for the lysis simulation.

This module defines all constants, enumerations, and default values used across
the lysis simulation package. It includes:

- Physical units management via Pint
- File naming conventions for input/output data
- Enumeration types for boundary conditions, fiber directions, molecule states, etc.
- Rectilinear grid neighbor calculation constants
- Random number draw type identifiers for Fortran compatibility
- NumPy dtype to format string mappings for text file I/O

The module provides a single ``CONST`` object that aggregates commonly used
constants for convenient access throughout the codebase.
"""

from enum import Enum, IntEnum, unique, Flag, auto

import numpy as np

from pint import UnitRegistry


# A Pint unit registry for use across the entire project.
# Any additional units should be added here.
# This registry handles unit conversions and dimensional analysis.
ureg = UnitRegistry()
# Shorthand for creating Quantity objects with units
Q_ = ureg.Quantity


class Const:
    """Aggregator class for commonly used constants and enumerations.

    Provides a single object that collects all major enumerations and constant
    collections for convenient access throughout the codebase. This eliminates
    the need to import multiple enum classes separately.

    :ivar BOUND: Enumeration of boundary directions (TOP, BOTTOM, LEFT, RIGHT, etc.)
    :vartype BOUND: BoundaryDirection
    :ivar BOUND_COND: Enumeration of boundary condition types (REFLECTING, PERIODIC, etc.)
    :vartype BOUND_COND: BoundaryCondition
    :ivar DIR: Enumeration of fiber directions in the hexagonal grid
    :vartype DIR: FiberDirection
    :ivar NEIGHBORHOOD: Constants for calculating neighbors in the edge grid
    :vartype NEIGHBORHOOD: Neighbors
    :ivar MOL_STATUS: Enumeration of molecule binding states
    :vartype MOL_STATUS: MolStatus
    :ivar DATASET_STORAGE_TYPE: Enumeration of data storage formats
    :vartype DATASET_STORAGE_TYPE: DataSetStorageType
    :ivar TPA_LEAVE_TIME_BINS: Number of tPA-leaving-time quantile bins the
        micro->macro conversion builds degradation-time distributions for (100)
    :vartype TPA_LEAVE_TIME_BINS: int
    :ivar NUMPY_SAVETXT_FORMATS: Format strings for numpy.savetxt by dtype
    :vartype NUMPY_SAVETXT_FORMATS: dict

    Example:
        >>> from lysis.config.constants import CONST
        >>> if boundary == CONST.BOUND.TOP:
        >>>     # Handle top boundary
        >>>     pass
    """

    def __init__(self):
        self.BOUND = BoundaryDirection
        self.BOUND_COND = BoundaryCondition
        self.DIR = FiberDirection
        self.NEIGHBORHOOD = Neighbors()
        self.MOL_STATUS = MolStatus
        self.DATASET_STORAGE_TYPE = DataSetStorageType
        # Number of bins the micro->macro conversion carves the tPA
        # leaving-time distribution into.  generate_macroscale_in() sorts the
        # microscale simulations by tPA leaving time and splits them into this
        # many equal-count quantile bins; each bin yields one sorted
        # distribution of fiber degradation times (one column of lysismat /
        # binned_fiber_degrade_time).  The macroscale model later draws a
        # uniform r in [0, 1] and multiplies by this value to pick which bin's
        # degradation distribution to sample.  Consequently:
        #   * lysismat / binned_fiber_degrade_time have this many columns,
        #   * lenlysisvect / binned_fiber_degraded have this many entries,
        #   * tPAleave / tsectPA / bin_edge_* have this many + 1 (the bin edges),
        #   * the Fortran macro binary's --nummicro is micro_simulations // this.
        self.TPA_LEAVE_TIME_BINS = 100
        self.DATASPEC_VERSION_ATTR = "dataspec_version"
        self.CONVERTED_FROM_ATTR = "converted_from"
        self.RENAMED_FROM_ATTR = "renamed_from"
        # Init-time provenance (stamped by init-experiment / init-macroscale).
        # Records the ``src/lysis/`` Python layer at init time.  ``init_dirty``
        # is the 3-state string ``"clean"`` / ``"dirty"`` / ``"unknown"``.
        self.INIT_VERSION_ATTR = "init_version"
        self.INIT_DIRTY_ATTR = "init_dirty"
        self.INIT_TIMESTAMP_ATTR = "init_timestamp"
        self.INIT_HOSTNAME_ATTR = "init_hostname"
        # Pipeline-time provenance (stamped by run-micro / run-macro at HDF5
        # import).  Records the ``src/lysis/`` Python pipeline that orchestrated
        # init/import -- NOT the engine that executes the simulation (that is
        # the backend, below).  ``pipeline_dirty`` is the 3-state string
        # ``"clean"`` / ``"dirty"`` / ``"unknown"``.
        self.PIPELINE_VERSION_ATTR = "pipeline_version"
        self.PIPELINE_DIRTY_ATTR = "pipeline_dirty"
        self.PIPELINE_TIMESTAMP_ATTR = "pipeline_timestamp"
        self.PIPELINE_HOSTNAME_ATTR = "pipeline_hostname"
        # Backend provenance -- the simulation engine that actually produced
        # the output.  For Fortran runs these come from the binary's
        # ``--version`` output; ``backend_commit`` / ``backend_dirty`` /
        # ``backend_compiler`` / ``backend_type`` are stamped unconditionally
        # by ``run-*`` commands.  ``backend_dirty`` is the 3-state string
        # ``"clean"`` / ``"dirty"`` / ``"unknown"``.  ``backend_type`` is
        # ``"fortran"`` (or ``"python"`` for a future Python backend).
        # ``stale_backend_override`` is set only when the binary's stamp
        # disagrees with the source tree and the user has explicitly
        # overridden the preflight check.
        self.BACKEND_COMMIT_ATTR = "backend_commit"
        self.BACKEND_DIRTY_ATTR = "backend_dirty"
        self.BACKEND_COMPILER_ATTR = "backend_compiler"
        self.BACKEND_TYPE_ATTR = "backend_type"
        self.STALE_BACKEND_OVERRIDE_ATTR = "stale_backend_override"
        # Boolean marker stamped (``True``) by the historical-build workflow
        # (``lysis run-{micro,macro} --fortran-commit <ref>``); absent on every
        # other run.  Lets a reader of the HDF5 file see at a glance that this
        # run's binary was rebuilt from an earlier commit rather than being the
        # build that matches the current source tree.  The commit SHA itself
        # lives in ``backend_commit``.
        self.BACKEND_HISTORICAL_ATTR = "backend_historical"
        # Records the RNG seeding scheme when it departs from the default
        # SeedSequence "split".  Stamped (``"direct"``) by ``run-micro
        # --direct`` for legacy reproduction, where the seed is fed straight
        # to the Fortran KISS RNG with no SeedSequence interposition.  Absent
        # ⇒ "split".  This is a provenance stamp, not a dataspec dataset.
        self.SEED_SCHEME_ATTR = "seed_scheme"
        self.LYSIS_ALLOW_STALE_BINARY_ENV = "LYSIS_ALLOW_STALE_BINARY"
        self.LYSIS_ALLOW_DIRTY_ENV = "LYSIS_ALLOW_DIRTY"
        self.LYSIS_ALLOW_COMMIT_MISMATCH_ENV = "LYSIS_ALLOW_COMMIT_MISMATCH"
        self.NUMPY_SAVETXT_FORMATS = {
            # Boolean types - save as 0 or 1
            np.dtype("bool"): "%d",
            # Signed integer types
            np.dtype("int8"): "%d",
            np.dtype("int16"): "%d",
            np.dtype("int32"): "%d",
            np.dtype("int64"): "%d",
            # Unsigned integer types
            np.dtype("uint8"): "%u",
            np.dtype("uint16"): "%u",
            np.dtype("uint32"): "%u",
            np.dtype("uint64"): "%u",
            # Floating point types - use scientific notation with full precision
            np.dtype("float32"): "%.9e",  # 32-bit float: ~7-9 decimal digits
            np.dtype("float64"): "%.18e",  # 64-bit float: ~15-18 decimal digits
            # Object and string types - use string representation
            np.dtype("object"): "%s",
            np.dtype("U0"): "%s",  # Unicode strings (any length)
        }

    def get_savetxt_format(self, dtype):
        """Get the numpy.savetxt format string for a given dtype.

        This method handles dtype matching by trying exact match first, then
        falling back to dtype kind matching for string types and structured arrays.

        :param dtype: NumPy dtype to get format for
        :type dtype: numpy.dtype or str
        :return: Format string for numpy.savetxt
        :rtype: str
        :raises ValueError: If dtype is not supported

        Example:
            >>> from lysis.config.constants import CONST
            >>> import numpy as np
            >>> fmt = CONST.get_savetxt_format(np.dtype('float64'))
            >>> print(fmt)  # '%.18e'
            >>> fmt = CONST.get_savetxt_format(np.dtype('<U75'))
            >>> print(fmt)  # '%s'
        """
        dtype = np.dtype(dtype)

        # Try exact match first
        if dtype in self.NUMPY_SAVETXT_FORMATS:
            return self.NUMPY_SAVETXT_FORMATS[dtype]

        # Handle Unicode strings of any length (e.g., '<U75')
        if dtype.kind == "U":
            return self.NUMPY_SAVETXT_FORMATS[np.dtype("U0")]

        # Handle byte strings (e.g., 'S10')
        if dtype.kind == "S":
            return "%s"

        # Handle structured arrays - return format for each field
        if dtype.names is not None:
            # Structured dtype - create format list for all fields
            formats = []
            for name in dtype.names:
                field_dtype = dtype.fields[name][0]
                formats.append(self.get_savetxt_format(field_dtype))
            return formats

        # Unknown dtype
        raise ValueError(
            f"No numpy.savetxt format defined for dtype: {dtype}\n"
            f"Consider adding it to CONST.NUMPY_SAVETXT_FORMATS"
        )


class Neighbors:
    """Constants for calculating neighboring edges in a rectilinear edge grid.

    This class provides the coordinate deltas needed to find the 8 neighbors of
    any edge in the rectilinear grid. The grid represents fibrin fibers as edges
    in a rectilinear lattice, where each edge has up to 8 neighboring edges.

    The grid uses three edge types (X, Y, Z) corresponding to the three edge
    orientations in a rectilinear lattice. Each edge type has different neighbor
    offset patterns.

    Boundary reflection constants handle the special cases when an edge is at
    the grid boundary and periodic/reflecting boundary conditions apply.

    :ivar X: Coordinate deltas for X-type edges: (row_deltas, col_deltas)
    :vartype X: Tuple[Tuple[int], Tuple[int]]
    :ivar Y: Coordinate deltas for Y-type edges: (row_deltas, col_deltas)
    :vartype Y: Tuple[Tuple[int], Tuple[int]]
    :ivar Z: Coordinate deltas for Z-type edges: (row_deltas, col_deltas)
    :vartype Z: Tuple[Tuple[int], Tuple[int]]
    :ivar TOP_REFL: Additional row deltas when at top boundary
    :vartype TOP_REFL: Tuple[int]
    :ivar BOTTOM_REFL: Additional row deltas when at bottom boundary
    :vartype BOTTOM_REFL: Tuple[int]
    :ivar LEFT_REFL: Additional column deltas when at left boundary
    :vartype LEFT_REFL: Tuple[int]
    :ivar RIGHT_REFL: Additional column deltas when at right boundary
    :vartype RIGHT_REFL: Tuple[int]

    Note:
        See EdgeGrid.neighbor() method documentation for details on how these
        constants are applied to calculate neighbor coordinates.
    """

    def __init__(self):
        # The changes that need to be made to the 2-D coordinates to obtain a neighbor.
        # For example, to obtain the first neighbor of an X-edge,
        # you decrease the row coordinate by 1 and decrease the column coordinate by 2.
        # Each tuple contains 8 values (one for each of the 8 possible neighbors)
        self.X = ((-1, -1, 0, 0, 0, 0, 0, 0), (-2, 1, -2, 1, -1, -1, 2, 2))
        self.Y = ((0, 0, 1, 1, 0, 1, 0, 1), (1, 1, 1, 1, -1, -1, 2, 2))
        self.Z = ((-1, -1, 0, 0, 0, 0, 0, 0), (-1, -1, -1, -1, -2, -2, 1, 1))
        # The additional changes that need to be made to the 2-D coordinates when at a boundary.
        # The top/bottom deltas need to be added to the row coordinate
        self.TOP_REFL = (0, 0, -1, -1, 0, 0, 0, 0)
        self.BOTTOM_REFL = (1, 1, 0, 0, 0, 0, 0, 0)
        # The left/right deltas need to be added to the column coordinate
        self.LEFT_REFL = (0, 0, 0, 0, 3, 3, 0, 0)
        self.RIGHT_REFL = (0, 0, 0, 0, 0, 0, -3, -3)


@unique
class RandomDraw(IntEnum):
    """Enumeration of random number draw types for Fortran compatibility.

    In Fortran duplication mode, the macroscale simulation pre-generates all
    random numbers for a timestep in a specific order to match the original
    Fortran implementation. This enum identifies which type of random draw
    each array position corresponds to.

    Each value represents a different stochastic event in the simulation that
    requires a random number draw.

    :cvar BINDING_TIME_WHEN_UNBINDING: Random number for calculating binding time after unbinding
    :cvar BINDING_TIME_WHEN_MOVING: Random number for calculating binding time after moving
    :cvar MICRO_UNBIND: Random number for determining if unbinding is forced (micro unbind)
    :cvar MOVE: Random number for movement decisions
    :cvar UNBINDING_TIME: Random number for selecting unbinding time from distribution
    :cvar LYSIS_TIME: Random number for selecting lysis time from distribution
    :cvar CONFLICT_RESOLUTION: Random number for resolving bind/move conflicts
    :cvar RESTRICTED_MOVE: Random number for restricted movement to degraded edges
    """

    BINDING_TIME_WHEN_UNBINDING = 0
    BINDING_TIME_WHEN_MOVING = 1
    MICRO_UNBIND = 2
    MOVE = 3
    UNBINDING_TIME = 4
    LYSIS_TIME = 5
    CONFLICT_RESOLUTION = 6
    RESTRICTED_MOVE = 7


@unique
class BoundaryDirection(IntEnum):
    """Enumeration of grid boundary directions.

    Identifies the six boundaries of the 3D simulation grid. Used to specify
    which boundary a boundary condition applies to or which boundary an edge/node
    is adjacent to.

    :cvar TOP: Top boundary (positive row direction)
    :cvar BOTTOM: Bottom boundary (negative row direction)
    :cvar LEFT: Left boundary (negative column direction)
    :cvar RIGHT: Right boundary (positive column direction)
    :cvar FRONT: Front boundary (negative depth direction, where tPA enters)
    :cvar BACK: Back boundary (positive depth direction, far side of clot)
    """

    TOP = 0
    BOTTOM = 1
    LEFT = 2
    RIGHT = 3
    FRONT = 4
    BACK = 5


@unique
class BoundaryCondition(Enum):
    """Enumeration of boundary condition types.

    Defines how molecules behave when they reach the edge of the simulation grid.
    Different boundary conditions can be applied to different boundaries.

    :cvar REFLECTING: Molecules bounce back when hitting the boundary (elastic reflection)
    :cvar PERIODIC: Molecules wrap around to the opposite boundary (toroidal topology)
    :cvar CONTINUING: Molecules can pass through the boundary freely (open boundary)
    """

    REFLECTING = 0
    PERIODIC = 1
    CONTINUING = 2


@unique
class FiberDirection(Enum):
    """Enumeration of fiber edge directions in the rectilinear grid.

    Identifies the orientation of fiber edges in the 3D rectilinear lattice.
    Each fiber edge runs in one of six directions, and opposite directions
    are encoded with opposite signs for convenient direction reversal.

    The three pairs of opposite directions correspond to the three types of
    edges (X, Y, Z) in the rectilinear grid structure.

    :cvar UP: Upward vertical direction (positive row, Z-type edge)
    :cvar DOWN: Downward vertical direction (negative row, Z-type edge)
    :cvar LEFT: Leftward horizontal direction (negative column, X-type edge)
    :cvar RIGHT: Rightward horizontal direction (positive column, X-type edge)
    :cvar OUT: Outward depth direction (positive depth, Y-type edge)
    :cvar IN: Inward depth direction (negative depth, Y-type edge)

    Note:
        The numeric values' signs indicate opposite directions, allowing
        direction reversal via negation: -FiberDirection.UP == FiberDirection.DOWN
    """

    UP = 1
    DOWN = -1
    LEFT = 2
    RIGHT = -2
    OUT = 3
    IN = -3


@unique
class RunComponent(Flag):
    """Flags for specifying which simulation components to execute.

    This Flag enum uses bit flags to allow combining multiple components using
    bitwise OR operations. This enables flexible specification of which parts
    of the simulation pipeline should run.

    The simulation pipeline consists of:
    1. Microscale simulation (generate binding/unbinding statistics)
    2. Microscale postprocessing (convert microscale output for macroscale input)
    3. Macroscale simulation (simulate tPA diffusion and fiber degradation)
    4. Macroscale postprocessing (analyze and visualize macroscale results)

    :cvar NONE: Run no components (placeholder)
    :cvar MICRO: Run microscale simulation
    :cvar MICRO_POSTPROCESSING: Run microscale postprocessing
    :cvar MACRO: Run macroscale simulation
    :cvar MACRO_POSTPROCESSING: Run macroscale postprocessing
    :cvar ALL: Run all components (bitwise OR of all components)

    Example:
        >>> # Run only microscale components
        >>> components = RunComponent.MICRO | RunComponent.MICRO_POSTPROCESSING
        >>> # Check if macro should run
        >>> if RunComponent.MACRO in components:
        >>>     run_macroscale()
    """

    NONE = 0
    MICRO = 1
    MICRO_POSTPROCESSING = 2
    MACRO = 4
    MACRO_POSTPROCESSING = 8
    ALL = MICRO | MICRO_POSTPROCESSING | MACRO | MACRO_POSTPROCESSING


@unique
class MolStatus(IntEnum):
    """Enumeration of tPA molecule binding states.

    Tracks the current state of tPA molecules in the macroscale simulation.
    Molecules transition between these states based on binding/unbinding events
    and fiber degradation.

    - **Bound** molecules are bound to an intact part of the fibrin lattice.
    - **Unbound** molecules are not bound to any fibrin, either because they
      have never bound, or because they kinetically unbound from the fibrin.
      This kinetic unbinding takes place at the microscale level.
    - **Micro-unbound** (forced unbinding) molecules are still bound to a
      small piece of fibrin, but that piece has been separated from the lattice
      due to plasmin-mediated degradation of the binding site they occupy. This
      occurs at the microscale level. The molecule cannot rebind while in this
      state, but the fibrin fragment it is attached to might move.
    - **Macro-unbound** (unbinding by degradation) molecules are still bound
      to a large piece of fibrin, but that piece has been separated from the
      lattice because the fiber has passed the ``snap_percentage`` threshold
      and is considered fully degraded. This occurs at the macroscale level.
      The molecule cannot rebind while in this state, but the fibrin fragment
      it is attached to might move.

    Most code follows the "into-and-along" schema, by which unbound molecules
    are free to move anywhere on the edge grid and are free to bind to any
    intact fiber. Micro-unbound molecules also have unrestricted movement,
    since the fibrin fragment they are attached to is small enough to pass
    through the pores in the lattice, but they cannot rebind until their
    waiting period has elapsed. During their waiting period, macro-unbound
    molecules have restricted movement and cannot rebind. Restricted movement
    means the molecule can only move to edges where no intact fiber exists
    (either in the fibrin-free region, or where the fiber has fully degraded).
    This is because the fibrin fragment they are attached to is too large to
    pass through intact parts of the lattice.

    For details, see
    Bannish, B. E., Paynter, B., Risman, R. A., Shroff, M., & Tutwiler, V. (2024).
    The effect of plasmin-mediated degradation on fibrinolysis and
    tissue plasminogen activator diffusion. Biophysical Journal, 123(5), 610-621.

    :cvar UNBOUND: Molecule is unbound and free to move and bind normally.
    :cvar BOUND: Molecule is bound to a fiber. It cannot move or bind to any
        other fiber.
    :cvar MACRO_UNBOUND: Unbound by fiber degradation. Restricted movement
        and cannot rebind until waiting period expires.
    :cvar MICRO_UNBOUND: Forced unbinding at the microscale. Unrestricted
        movement but cannot rebind until waiting period expires.
    """

    UNBOUND = 0
    BOUND = 1
    MACRO_UNBOUND = 2
    MICRO_UNBOUND = 3


@unique
class DataSetStorageType(Enum):
    """Enumeration of data storage formats and locations.

    Specifies how and where simulation data should be stored. The simulation
    supports multiple storage backends to accommodate different data sizes,
    access patterns, and interoperability requirements.

    HDF5-based storage options are preferred for large numerical datasets due
    to efficient compression and partial I/O. File-based options are useful for
    smaller datasets, human-readable output, or legacy format compatibility.

    :cvar HDF5_ATTR: Store as HDF5 attribute (small metadata, <64KB)
    :cvar HDF5_GROUP: Store as HDF5 group (organizational container)
    :cvar HDF5_DATASET: Store as HDF5 dataset (large numerical arrays, primary data storage)
    :cvar FILE_HDF5: Store as standalone HDF5 file
    :cvar FILE_TEXT: Store as human-readable text file
    :cvar FILE_PARSED: Store as human-readable text file that must be parsed for values
    :cvar FILE_BINARY: Store as binary file (compact but not human-readable)
    :cvar FILE_JSON: Store as JSON file (human-readable, good for metadata)

    Note:
        HDF5 options require the h5py library. File options use standard Python I/O.
    """

    HDF5_ATTR = auto()
    HDF5_GROUP = auto()
    HDF5_DATASET = auto()
    FILE_HDF5 = auto()
    FILE_TEXT = auto()
    FILE_PARSED = auto()
    FILE_BINARY = auto()
    FILE_JSON = auto()


#: Global constants object for convenient access to all enumerations.
#:
#: This singleton instance provides access to all constant enumerations through
#: a single import. Use this instead of importing individual enum classes.
#:
#: Example:
#:     >>> from lysis.config.constants import CONST
#:     >>> if molecule_state == CONST.MOL_STATUS.BOUND:
#:     >>>     process_bound_molecule()
#: Mapping of fiber type codes to their physical parameters.
#:
#: These codes identify standard fiber configurations used in the fibrinolysis
#: simulation.  Each entry maps a short code (from the header of
#: ``micro_rates.f90``) to the fiber radius and the number of protofibril
#: nodes per row in the microscale lattice.
#:
#: The fiber radius is half the fiber bundle diameter listed in the source
#: code.  The nodes-per-row value determines the ``nodes_in_micro_row``
#: parameter in :class:`~lysis.config.parameters.MicroParameters`.
#:
#: Example:
#:     >>> from lysis.config.constants import FIBER_TYPES
#:     >>> FIBER_TYPES["Q4"]["nodes_in_micro_row"]
#:     13
FIBER_TYPES = {
    "Q0":     {"fiber_radius": Q_("23.0 nanometers"),  "nodes_in_micro_row": 4},
    "Q1":     {"fiber_radius": Q_("28.7 nanometers"),  "nodes_in_micro_row": 5},
    "Q2":     {"fiber_radius": Q_("36.35 nanometers"), "nodes_in_micro_row": 7},
    "Q3":     {"fiber_radius": Q_("40.65 nanometers"), "nodes_in_micro_row": 8},
    "TF-v":   {"fiber_radius": Q_("52.55 nanometers"), "nodes_in_micro_row": 5},
    "TF-vii": {"fiber_radius": Q_("52.55 nanometers"), "nodes_in_micro_row": 7},
    "TF-x":   {"fiber_radius": Q_("52.55 nanometers"), "nodes_in_micro_row": 10},
    "TB-xi":  {"fiber_radius": Q_("61.5 nanometers"),  "nodes_in_micro_row": 11},
    "TB-xiii": {"fiber_radius": Q_("61.5 nanometers"), "nodes_in_micro_row": 13},
    "Q4":     {"fiber_radius": Q_("72.7 nanometers"),  "nodes_in_micro_row": 13},
}

#: Fortran sentinel (v1.90.0 and later) marking a fibrin edge on which no tPA
#: molecule has yet landed, so no degradation has been scheduled for it.
#:
#: Set by ``t_degrade = 9.9d+100`` in
#: ``src/fortran/macro_diffuse_into_and_along__internal.f90``.
UNSCHEDULED_DEGRADE_TIME: float = 9.9e100

#: Degrade time of an empty (ghost) edge -- one of the ``empty_edges`` edges
#: below the clot, which can never degrade.
#:
#: Set by ``t_degrade(1:enoFB) = 0`` in
#: ``src/fortran/macro_diffuse_into_and_along__internal.f90``.
EMPTY_EDGE_DEGRADE_TIME: float = 0.0

#: The v1.85.0 "no degradation scheduled" sentinel.
#:
#: v1.85.0 initialised the *whole* ``t_degrade`` vector to zero
#: (``t_degrade=0.0d+00`` in ``macro_Q2_diffuse_along.f90``) and never
#: distinguished empty edges, so this value is ambiguous on its own: it means
#: "empty edge" for the first ``empty_edges`` entries and "fibrin edge, not yet
#: scheduled" for the rest.  Converters must disambiguate by index, never by
#: value.
V185_UNSCHEDULED_DEGRADE_TIME: float = 0.0

CONST = Const()
