"""Rectilinear edge grid representation for fibrin network simulation.

This module provides the EdgeGrid class and associated utilities for representing
and manipulating a 3D rectilinear lattice of fibrin fiber edges. The grid is the
fundamental spatial structure for the macroscale simulation of tPA-mediated clot
lysis.

Grid Structure
--------------

The EdgeGrid represents fibrin fibers as edges in a rectilinear lattice. Each edge
connects two nodes and represents a single fiber segment. The grid is organized into
rows (y-direction), with each row containing nodes in the x-direction, and depth
in the z-direction.

Three edge types exist, corresponding to the three spatial directions:

- **x-edges**: Horizontal edges (left-right), rank j where j % 3 == 2
- **y-edges**: Vertical edges (up-down), rank j where j % 3 == 0
- **z-edges**: Depth edges (in-out), rank j where j % 3 == 1

Each edge has 8 neighbors, pre-computed for efficient access during simulation.

Indexing Systems
----------------

The module supports two indexing systems:

1. **Python 2D indexing (i, j)**:
   - i: row index (0 to rows-1)
   - j: rank within row (0 to 3*cols-2)
   - Used throughout Python implementation
   - Edges ordered as triplets: (y-edge, z-edge, x-edge) for each node

2. **Fortran 1D indexing**:
   - Single index from 0 to total_edges-1
   - Orders edges differently: all z/x edges in a row, then all y edges
   - Used for compatibility with legacy Fortran code
   - **Fortran uses 1-based indexing; add 1 to all indices for Fortran**

Conversion functions (from_fortran_edge_index, to_fortran_edge_index) handle
translation between these systems.

Boundary Conditions
-------------------

The grid supports three boundary condition types (from constants.BoundaryCondition):

- **REFLECTING**: Molecules bounce back at boundaries
- **PERIODIC**: Molecules wrap around to opposite boundary
- **CONTINUING**: Molecules pass through freely (open boundary)

The top row has no y-edges due to the reflecting top boundary condition. Empty
rows at the front represent space where tPA starts before entering the clot.

Main Classes
------------

- **EdgeGrid**: Main class representing the grid, managing fiber status, and
  computing neighbors with boundary conditions

Utility Functions
-----------------

- **from_fortran_edge_index**: Convert Fortran 1D index to Python (i,j)
- **to_fortran_edge_index**: Convert Python (i,j) to Fortran 1D index
- **from_fortran_edge_index_array**: Vectorized Fortran to Python conversion
- **to_fortran_edge_index_array**: Vectorized Python to Fortran conversion
- **generate_fortran_neighborhood_structure**: Generate neighbor lookup table for Fortran

Example Usage
-------------

Create an edge grid for a simulation::

    >>> from lysis.util import Run
    >>> run = Run(...)  # Configure simulation parameters
    >>> grid = EdgeGrid(run)
    >>>
    >>> # Get the 3rd neighbor of edge at (10, 5)
    >>> neighbor_i, neighbor_j = grid.neighbor(10, 5, 3)
    >>>
    >>> # Check fiber status
    >>> if grid.fiber_status[10, 5] > current_time:
    >>>     print("Fiber is intact")

Pre-compute all neighbors for fast lookup::

    >>> neighbors = EdgeGrid.generate_neighborhood_structure(run)
    >>> # neighbors[edge_idx, k] gives the 1D index of the k-th neighbor

Convert between indexing systems::

    >>> fortran_idx = to_fortran_edge_index(10, 5, rows=100, nodes_in_row=50)
    >>> i, j = from_fortran_edge_index(fortran_idx, rows=100, nodes_in_row=50)

See Also
--------

- constants.py: Neighbor offset constants and boundary condition enumerations
- np_macroscale.py: Macroscale simulation that uses EdgeGrid
"""

from functools import partial
from typing import Tuple

import numpy as np

from pint import Quantity

from ..config.constants import Const, BoundaryCondition
from ..execution.run import Run


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

CONST = Const()


class EdgeGrid(object):
    """The main class containing a 3-D grid of edges. This represents an
    xy-planar slice, one edge high, of a clot.
    An edge's location is determined by its Row and its Rank within that row.

    Co-ordinate arrangement::

        '         {        /           /           /           /           /
        '   Row 3 {       1           4           7           10          13
        '         {      /           /           /           /           /
        '         {     +-----2-----+-----5-----+-----8-----+-----11----+
        '              /|          /|          /|          /|          /|
        '             / |         / |         / |         / |         / |
        '         {     0  /        3  /        6  /        9  /       12  /
        '   Row 2 {     | 1         | 4         | 7         | 10        | 13
        '         {     |/          |/          |/          |/          |/
        '         {     +-----2-----+-----5-----+-----8-----+-----11----+
        '              /|          /|          /|          /|          /|
        '             / |         / |         / |         / |         / |
        '         {     0  /        3  /        6  /        9  /       12  /
        '   Row 1 {     | 1         | 4         | 7         | 10        | 13
        '         {     |/          |/          |/          |/          |/
        '         {     +-----2-----+-----5-----+-----8-----+-----11----+
        '              /|          /|          /|          /|          /|
        '             / |         / |         / |         / |         / |
        '         {     0  /        3  /        6  /        9  /       12  /
        '   Row 0 {     | 1         | 4         | 7         | 10        | 13
        '         {     |/          |/          |/          |/          |/
        '         {     +-----2-----+-----5-----+-----8-----+-----11----+
        '              /           /           /           /           /
        '             /           /           /           /           /

    """

    def __init__(
        self,
        run: Run,
        boundary_conditions: Tuple[BoundaryCondition, BoundaryCondition] | None = None,
        initial_fiber_status: float = float("inf"),
    ):
        """Initialize an EdgeGrid for the rectilinear lattice simulation.

        :param run: The Run object containing simulation parameters. This provides
            the grid dimensions (rows, cols) and empty row configuration.
        :type run: Run
        :param boundary_conditions: Tuple of (top/bottom, left/right) boundary
            conditions. If None, defaults to reflecting boundaries on all sides.
        :type boundary_conditions: Tuple[BoundaryCondition, BoundaryCondition], optional
        :param initial_fiber_status: Initial degradation time for all fibers.
            Should be infinity (or greater than simulation length) to indicate
            intact fibers. Empty rows are automatically set to 0.
        :type initial_fiber_status: float, optional
        """

        self.total_rows = run.macro_params.rows
        """int: The total number of rows in the grid."""
        self.nodes_in_row = run.macro_params.cols
        """int: The number of nodes (not edges) in each row of this EdgeGrid."""
        self.empty_rows = run.macro_params.empty_rows
        """int: The number of empty (fibrin-free) rows in the grid."""
        # Set the appropriate boundary conditions
        if boundary_conditions is not None:
            self.boundary_conditions = boundary_conditions
        else:
            self.boundary_conditions = (
                CONST.BOUND_COND.REFLECTING,
                CONST.BOUND_COND.REFLECTING,
            )

        self.ranks = 3 * self.nodes_in_row - 1
        """int: The number of fibers in one row of this grid."""
        self.fiber_rows = self.total_rows - self.empty_rows
        """int: The nunber of rows of the grid that contain actual fibers."""
        self.fiber_status = initial_fiber_status * np.ones(
            (self.total_rows, self.ranks), dtype=np.double
        )
        """np.ndarray: The status of the fibers in this EdgeGrid. This is 
                    essentially the degrade time of the fiber.
                    If degrade time < current time, the fiber is degraded.
                    If degrade time = 0, the edge contains no fiber
                    e.g., self.fiber_status[60, 25] is the status of the fiber 
                    located along edge [60, 25].

                    Supports two-dimensional slicing."""
        # Set the status of empty rows to zero. 
        # That is, the "fiber" here "degraded" at time zero 
        # (since it never actually existed in the first place).
        self.fiber_status[: self.empty_rows] = 0

    def _is_valid_index(self, i: int, j: int) -> str | None:
        """Check if [i, j] is a valid index for this EdgeGrid.

        Validates that the given coordinates refer to an existing edge in the grid,
        accounting for grid boundaries and the fact that y-edges don't exist on
        the top row with reflecting boundary conditions.

        :param i: The row index of the edge (0 to total_rows-1)
        :type i: int
        :param j: The rank (column) index of the edge within its row (0 to ranks-1)
        :type j: int
        :return: None if the index is valid, otherwise an error message string
            describing why the index is invalid
        :rtype: str or None
        """
        # Edges are 0 through self.ranks-1
        if j < 0 or j > self.ranks - 1:
            return (
                f"Index j={j} out of bounds. "
                f"This model only has edges [0..{self.ranks - 1}] in "
                f"each row."
            )
        # Rows are 0 through self.rows-1
        if i < 0 or i > self.total_rows - 1:
            return (
                f"Index i={i} out of bounds. This model only has rows "
                f"[0..{self.total_rows - 1}]."
            )
        # If this EdgeGrid is at the top of a larger sliced grid, or is the
        # whole grid itself, then the top row has no y-edges in it.
        if (
            i == self.total_rows - 1
            and self.boundary_conditions[CONST.BOUND.TOP] == CONST.BOUND_COND.REFLECTING
            and j % 3 == 0
        ):
            return (
                f"y-edges do not exist on the top row of this grid "
                f"(row {self.total_rows - 1}). Location ({i}, {j})"
            )
        # Everything seems fine, so return None
        return None

    def neighbor(self, i: int, j: int, k: int) -> Tuple[int, int]:
        """Find the coordinates of a neighboring edge.

        Given the coordinates of an edge and the index of a neighbor in its
        neighborhood, returns the coordinates of that neighbor, accounting for
        boundary conditions.

        Each edge has 8 neighbors arranged according to its type (x, y, or z edge).
        The edge type is determined by j % 3: 0=y-edge, 1=z-edge, 2=x-edge.

        Neighborhoods are labeled with k:(i+p,j+q) where (p,q) is the offset
        from the generating edge to the k-th neighbor::

            * x-edge neighborhood::

                '        2:(i, j-2)      3:(i, j+1)
                '            |              |
                '            |  /           |  /
                '            | 5:(i, j-1)   | 7:(i, j+2)
                '            |/             |/
                '            +----(i, j)----+
                '          / |             /|
                ' 4:(i, j-1) |   6:(i, j+2) |
                '        /   |          /   |
                '            |              |
                '      0:(i-1, j-2)   1:(i-1, j+1)

            * y-edge neighborhood::

                '                     /
                '                    3:(i+1, j+1)
                '                   /
                ' 5:(i+1, j-1)-----+-----7:(i+1, j+2)
                '                 /|
                '     2:(i+1, j+1) |
                '               /  |
                '               (i, j)
                '                  |  /
                '                  | 1:(i, j+1)
                '                  |/
                '   4:(i, j-1)-----+-----6:(i, j+2)
                '                 /
                '       0:(i, j+1)
                '               /

            * z-edge neighborhood::

                '                          |
                '               3:(i+1, j-1)
                '                          |
                '         5:(i+1, j-2)-----+-----7:(i+1, j+1)
                '                         /|
                '                        / 1:(i, j-1)
                '                       /  |
                '                   (i, j)
                '                  |  /
                '       2:(i+1, j-1) /
                '                  |/
                ' 4:(i+1, j-2)-----+-----6:(i+1, j+1)
                '                  |
                '                  0:(i, j-1)
                '                  |

        :param i: The row index of the generating edge (0 to total_rows-1)
        :type i: int
        :param j: The rank index of the generating edge within its row (0 to ranks-1)
        :type j: int
        :param k: The neighbor index in the neighborhood (0 to 7)
        :type k: int
        :return: The (row, rank) coordinates of the k-th neighbor
        :rtype: Tuple[int, int]
        :raises IndexError: If i, j, or k are out of bounds, or if (i,j) refers
            to a non-existent edge
        """
        # Check that we are in-bounds

        # Neighbors are 0 through 7
        if k < 0 or k > 7:
            raise IndexError(
                f"Index k={k} out of bounds. Each neighborhood has items [0..7]."
            )
        # Check if i and j are valid
        valid_index = self._is_valid_index(i, j)
        if valid_index is not None:
            raise IndexError(valid_index)

        # The index of the neighboring fiber being requested
        neighbor_i = int(i)
        neighbor_j = int(j)

        # Move to the required neighbor
        if j % 3 == 0:  # We are a y-edge
            neighbor_i += CONST.NEIGHBORHOOD.Y[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.Y[1][k]
        elif j % 3 == 1:  # We are a z-edge
            neighbor_i += CONST.NEIGHBORHOOD.Z[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.Z[1][k]
        elif j % 3 == 2:  # We are an x-edge
            neighbor_i += CONST.NEIGHBORHOOD.X[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.X[1][k]

        # Deal with boundary conditions

        # The bottom boundary of the grid.
        # Note that, if the edge generating the neighborhood is a y-edge,
        # then its neighborhood only involves fibers on its own row, or the row
        # above. Thus, the neighborhood of y-edges never overruns the bottom of
        # the grid
        if (
            i == 0  # We are at the bottom of the grid
            and j % 3 > 0  # and it is a z- or x-edge
        ):
            if (
                self.boundary_conditions[CONST.BOUND.BOTTOM]
                == CONST.BOUND_COND.REFLECTING
            ):
                # If the bottom is reflective, then we simply shift any
                # neighbors on row -1, to the row above.
                neighbor_i += CONST.NEIGHBORHOOD.BOTTOM_REFL[k]

        # The top boundary of the grid.
        # Note that, if the edge generating the neighborhood is a y-edge,
        # then its neighborhood only involves fibers on its own row, or the row
        # above. But, if this row is the top of the entire simulation
        # (REFLECTING) then there are no y-edges on this row. If this row is
        # the top of one slice (CONTINUING) then this row represents the bottom
        # row of the next slice and should not be processed here
        elif (
            i == self.total_rows - 1  # We are at the top of the grid
            and j % 3 > 0  # and it is a z- or x-edge
        ):
            if self.boundary_conditions[CONST.BOUND.TOP] == CONST.BOUND_COND.REFLECTING:
                neighbor_i += CONST.NEIGHBORHOOD.TOP_REFL[k]

        # The left boundary of the grid.
        # Note that, if the edge generating the neighborhood is an x-edge,
        # Then its neighborhood never overruns the side of the grid
        if j <= 1:  # We are the left-most y- or z-edge
            neighbor_j += CONST.NEIGHBORHOOD.LEFT_REFL[k]

        # The right boundary of the grid.
        # Note that, if the edge generating the neighborhood is an x-edge,
        # Then its neighborhood never overruns the side of the grid
        elif j >= self.ranks - 2:  # We are the right-most y- or z-edge
            neighbor_j += CONST.NEIGHBORHOOD.RIGHT_REFL[k]

        # Return the co-ordinates of the requested neighbor.
        return np.uint32(neighbor_i), np.uint32(neighbor_j)

    @staticmethod
    def generate_neighborhood_structure(run: Run) -> np.ndarray[np.ushort]:
        """Generate a pre-computed neighbor lookup table for all edges.

        Creates a 2D array where each row contains the 1D indices of the 8 neighbors
        for the corresponding edge. This pre-computation significantly accelerates
        the macroscale simulation by eliminating repeated neighbor calculations.

        The 1D indexing used here is the flattened Python 2D index
        (via np.ravel_multi_index), NOT the Fortran 1D index. For Fortran
        compatibility, use generate_fortran_neighborhood_structure() instead.

        Boundary conditions are respected: edges at boundaries may have duplicate
        neighbors in their list (e.g., reflecting boundaries cause the edge to
        neighbor itself).

        :param run: The Run object containing grid parameters (rows, cols, boundary conditions)
        :type run: Run
        :return: Array of shape (total_edges, 8) where result[edge_idx, k] is the
            1D index of the k-th neighbor of edge_idx. Duplicates may appear due
            to boundary conditions.
        :rtype: np.ndarray[np.ushort]
        """
        # We humans would rather interact with the edge grid in 2-D co-ordinates, 
        # but it is faster to work with a 1-D array, 
        # so we use the NumPy function ravel_multi_index to do that for us.
        # We use functools.partial to pre-fill the grid dimensions
        edge_lookup = partial(
            np.ravel_multi_index,
            dims=(run.macro_params.rows, run.macro_params.full_row),
        )
        # Initialize the EdgeGrid for the given parameters
        edge_grid = EdgeGrid(run)
        # Initialize an empty array to contain the neighbor indices
        neighbors = np.empty(
            (run.macro_params.rows * run.macro_params.full_row, 8), dtype=np.ushort
        )
        # Loop over all edges in the grid and their eight neighbors.
        for i, j, k in np.ndindex(edge_grid.total_rows, edge_grid.ranks, 8):
            # If the edge is one of the non-existant vertical edges on the top row
            if i == edge_grid.total_rows - 1 and j % 3 == 0:
                # Set its neighbors to itself
                neighbors[
                    edge_lookup((i, j)),
                    k,
                ] = edge_lookup((i, j))
            else:
                # Else determine the neighbor from the EdgeGrid and convert it's index.
                neighbors[
                    edge_lookup((i, j)),
                    k,
                ] = edge_lookup(tuple(edge_grid.neighbor(i, j, k)))
        # Return the completed array
        return neighbors

    @staticmethod
    def generate_fortran_neighborhood_structure(run: Run) -> np.ndarray[int]:
        """Generate neighbor lookup table compatible with Fortran Macro code.

        Creates a neighbor structure using Fortran's 1D indexing and ordering
        conventions. This allows the Python code to generate input files for
        the legacy Fortran implementation.

        **Critical differences from generate_neighborhood_structure():**

        - **Row ordering**: Edges ordered by Fortran 1D index (z/x edges first, then y edges per row)
        - **Neighbor ordering**: Each row's neighbors sorted in ascending order (Fortran convention)
        - **Indexing**: Returns 0-indexed values; **must add 1 before passing to Fortran**

        The Python generate_neighborhood_structure() uses different ordering for both
        rows and neighbors, making the two arrays incompatible despite containing the
        same topological information.

        :param run: The Run object containing grid parameters
        :type run: Run
        :return: Array of shape (total_edges, 8) where result[fortran_edge_idx, k] is
            the 0-indexed Fortran 1D index of the k-th neighbor (in sorted order).
            **Add 1 to all values before using in Fortran.**
        :rtype: np.ndarray[int]

        Warning:
            The returned indices are 0-based. Fortran uses 1-based indexing, so add 1
            to all values before writing to files for Fortran consumption.
        """
        # Create the EdgeGrid object for this run
        edge_grid = EdgeGrid(run)
        # Initialize an empty array
        fort_neighbors = np.empty((run.macro_params.total_edges, 8), dtype="int")
        # Loop over all of the edges in the grid
        for f in range(run.macro_params.total_edges):
            # Identify the co-ordinates of the 2-D grid corresponding to this Fortran index
            i, j = from_fortran_edge_index(
                f, run.macro_params.rows, run.macro_params.cols
            )
            # Loop over the edge's neighbors
            for k in range(8):
                neighbor_i, neighbor_j = edge_grid.neighbor(i, j, k)
                # Convert the neighbor's 2-D index back to the 1-D Fortran index
                fort_neighbors[f, k] = to_fortran_edge_index(
                    neighbor_i, neighbor_j, run.macro_params.rows, run.macro_params.cols
                )
        # Sort the idices by row (which is how they appear in Fortran) and return
        return np.sort(fort_neighbors)

    @staticmethod
    def get_spatial_coordinates(i: int, j: int) -> Tuple[float, float, float]:
        """Calculate the spatial coordinates of an edge's center point.

        Returns the (x, y, z) coordinates of the center of an edge in the
        rectilinear node grid underlying the EdgeGrid. Edges connect nodes, so
        their centers are offset by 0.5 in the direction of the edge.

        Edge types (determined by j % 3):
        - j % 3 == 0: y-edge (vertical), centered at (x, y+0.5, z)
        - j % 3 == 1: z-edge (depth), centered at (x, y, z+0.5)
        - j % 3 == 2: x-edge (horizontal), centered at (x+0.5, y, z)

        :param i: The row index of the edge in the EdgeGrid
        :type i: int
        :param j: The rank index of the edge within its row
        :type j: int
        :return: The (x, y, z) coordinates of the edge's center in the node grid.
            Coordinates are in node-spacing units (not physical units).
        :rtype: Tuple[float, float, float]

        Todo:
            Adapt for non-square lattices (currently assumes cubic unit cells)
        """
        # TODO: Adapt for non-square lattices
        # Identify which node forms the first end of the edge
        x = j // 3
        y = i
        z = 0
        # Move half a step in the direction of the edge.
        match j % 3:
            case 0:
                y += 0.5
            case 1:
                z += 0.5
            case 2:
                x += 0.5
        # Return the tuple of co-ordinates
        return x, y, z

    @staticmethod
    def get_distance(
        run: Run, a: Tuple[int, int], b: Tuple[int, int], metric: str = "euclidian"
    ) -> Quantity:
        """
        Calculate the distance between the centers of two edges on an EdgeGrid. 
        The distance returned will be a Pint Quantity which includes appropriate units.

        :param run: The Run that this edge grid will be used for.
        :type run: Run
        :param a: The 2-D index of the first edge
        :type a: Tuple[int, int]
        :param b: The 2-D index of the second edge
        :type b: Tuple[int, int]
        :param metric: The metric being used to calculate the distance, defaults to "euclidian".
            Options are:
                - "euclidian": :math:`\sqrt{(x_1-x_2)^2+(y_1-y_2)^2+(z_1-z_2)^2}`
                - "manhattan": :math:`\lvert x_1-x_2\rvert+\lvert y_1-y_2\rvert+\lvert z_1-z_2\rvert`
                - "taxicab": :math:`\lvert x_1-x_2\rvert+\lvert y_1-y_2\rvert+\lvert z_1-z_2\rvert`
                - "2d_euclidian": :math:`\sqrt{(x_1-x_2)^2+(y_1-y_2)^2}`
        :type metric: str, optional
        :raises AttributeError: If the metric passed is not among those implemented.
        :return: A Pint Quantity object containing the distance between the centers of the two edges.
        :rtype: Quantity
        """
        # Get the location of the centers of the two edges
        a_coord = EdgeGrid.get_spatial_coordinates(*a)
        b_coord = EdgeGrid.get_spatial_coordinates(*b)
        if metric == "euclidian":
            # Sum the square differences of all three co-ordinates
            squares = sum((a_coord[k] - b_coord[k]) ** 2 for k in range(3))
            # Square-root the sum, multiply by the pore size, and return
            return run.macro_params.pore_size * squares**0.5
        elif metric in ["manhattan", "taxicab"]:
            # Sum the absolute differences of all three co-ordinates
            sides = sum(abs(a_coord[k] - b_coord[k]) for k in range(3))
            # Multiply by the pore size and return
            return run.macro_params.pore_size * sides
        if metric == "2d_euclidian":
            # Sum the square differences of just the x and y co-ordinates
            squares = sum((a_coord[k] - b_coord[k]) ** 2 for k in range(2))
            # Square-root the sum, multiply by the pore size, and return
            return run.macro_params.pore_size * squares**0.5
        else:
            raise AttributeError(f"{metric} metric not implemented yet.")

    @staticmethod
    def full_row(rows: int, nodes_in_row: int) -> int:
        """Calculate the number of edges in a full row of the edge grid.

        A full row (not the top row) contains 3 edges per node, except the last
        node which has no x-edge. This gives: 3 * nodes_in_row - 1 edges.

        The top row has no y-edges and should use xz_row() instead.

        :param rows: The number of rows in the edge grid (unused but kept for API consistency)
        :type rows: int
        :param nodes_in_row: The number of nodes in each row of the edge grid
        :type nodes_in_row: int
        :return: The number of edges in a full row: 3 * nodes_in_row - 1
        :rtype: int
        """
        # Three edges for each node except the last one, which has no x-edge
        return 3 * nodes_in_row - 1

    @staticmethod
    def xz_row(rows: int, nodes_in_row: int) -> int:
        """Calculate the number of x- and z-edges in a row (excluding y-edges).

        Each node contributes one x-edge and one z-edge, except the last node
        which has no x-edge. This gives: 2 * nodes_in_row - 1 edges.

        This is the edge count for the top row, which has no y-edges due to
        the reflecting boundary condition.

        :param rows: The number of rows in the edge grid (unused but kept for API consistency)
        :type rows: int
        :param nodes_in_row: The number of nodes in each row of the edge grid
        :type nodes_in_row: int
        :return: The number of x- and z-edges: 2 * nodes_in_row - 1
        :rtype: int
        """
        # Two edges for each node except the last one, which has no x-edge
        return 2 * nodes_in_row - 1

    @staticmethod
    def total_edges(rows: int, nodes_in_row: int) -> int:
        """Calculate the total number of edges in the entire edge grid.

        The total is computed as:
        - (rows - 1) full rows, each with full_row(rows, nodes_in_row) edges
        - Plus 1 top row with only x- and z-edges: xz_row(rows, nodes_in_row)

        The top row has no y-edges due to the reflecting boundary condition
        at the top of the grid.

        :param rows: The number of rows in the edge grid
        :type rows: int
        :param nodes_in_row: The number of nodes in each row of the edge grid
        :type nodes_in_row: int
        :return: The total number of edges in the grid
        :rtype: int
        """
        # Each row has a full set of edges, except the top row which has no y-edges
        return EdgeGrid.full_row(rows, nodes_in_row) * (rows - 1) + EdgeGrid.xz_row(
            rows, nodes_in_row
        )


def generate_fortran_neighborhood_structure(rows: int, nodes_in_row: int) -> np.ndarray:
    """Generate Fortran-compatible neighbor structure using vectorized operations.

    This is an optimized, standalone version of EdgeGrid.generate_fortran_neighborhood_structure()
    that uses vectorized NumPy operations for better performance. It produces the
    same neighbor structure but doesn't require creating an EdgeGrid object.

    Returns a 2D array where each row corresponds to an edge (ordered by Fortran 1D
    index), and each row contains the 0-indexed Fortran 1D indices of that edge's
    8 neighbors, sorted in ascending order.

    **Important**: The returned indices are 0-based. Add 1 before using in Fortran code.

    :param rows: The number of rows in the edge grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the edge grid
    :type nodes_in_row: int
    :return: Array of shape (total_edges, 8) where result[fortran_edge_idx, :] contains
        the sorted 0-indexed Fortran 1D indices of the 8 neighbors. Add 1 for Fortran.
    :rtype: np.ndarray

    Note:
        This function is significantly faster than the EdgeGrid method version for
        large grids due to vectorization, but produces identical results.
    """
    # Generate a list of all fortran-index edges and get the equivalent 2-d indeces
    edges = from_fortran_edge_index_array(
        np.arange(EdgeGrid.total_edges(rows, nodes_in_row)), rows, nodes_in_row
    )
    # Figure out which edges are x-, y-, and z-edges
    edge_type = edges[:, 1] % 3
    # Top edges are x- and z-edges with the largest row index
    top = (edges[:, 0] == rows - 1) & (edge_type != 0)
    # Bottom edges are x- and z-edges with the smallest row index
    bottom = (edges[:, 0] == 0) & (edge_type != 0)
    # Left edges are those with rank 0 or 1
    left = edges[:, 1] <= 1
    # Right edges are those with the largest 2 ranks
    right = edges[:, 1] >= EdgeGrid.full_row(rows, nodes_in_row) - 2
    # Add an axis in the middle
    edges = edges.reshape(-1, 1, 2)
    # and then duplicate 8 times along that axis.
    # This is where the 8 neighbors will go
    neighbors = np.repeat(edges, 8, axis=1)

    # Create a matrix for the neighbor delta
    adds = np.empty(neighbors.shape, dtype=int)
    # Get the appropriate neighbor delta for each edge type (x, y, z) and put them in the matrix
    adds[edge_type == 0] = np.array(CONST.NEIGHBORHOOD.Y).T
    adds[edge_type == 1] = np.array(CONST.NEIGHBORHOOD.Z).T
    adds[edge_type == 2] = np.array(CONST.NEIGHBORHOOD.X).T
    # Deal with the boundaries
    adds[top, :, 0] += CONST.NEIGHBORHOOD.TOP_REFL
    adds[bottom, :, 0] += CONST.NEIGHBORHOOD.BOTTOM_REFL
    adds[left, :, 1] += CONST.NEIGHBORHOOD.LEFT_REFL
    adds[right, :, 1] += CONST.NEIGHBORHOOD.RIGHT_REFL
    # Add the delta to the current location to get the actual neighbor indices
    # TODO: This casting should be fine, but it wouldn't hurt to add code to check it
    neighbors = np.add(neighbors, adds, out=neighbors, casting="unsafe")

    # Rearrange the axes so that its (8, -1, 2)
    neighbors = np.moveaxis(neighbors, [0, 1, 2], [1, 0, 2])
    # Slice up this matrix by neigbor number, convert its 2-D indices to 1-D,
    # and then stack them back together.
    neighbor_1d = np.array(
        [to_fortran_edge_index_array(x, rows, nodes_in_row) for x in neighbors]
    )
    # Cast, Transpose, sort, and return
    return np.sort(neighbor_1d.T)


def from_fortran_edge_index(
    index: int, rows: int, nodes_in_row: int
) -> Tuple[int, int]:
    """Convert a 1-dimensional Fortran edge index to 2-dimensional (i, j) coordinates.

    Converts from the 1D index used in Fortran "Macro" data structures and files
    to the (i, j) index used by this Python package.

    .. note::
        **All indices are zero-indexed.** Add 1 before using in Fortran code.

    This function aids compatibility with legacy Fortran code during the transition
    phase and should be deprecated once Python implementation is complete.

    The Fortran 1D indexing orders edges as: all z-edges and x-edges in row 0,
    then all y-edges in row 0, then all z-edges and x-edges in row 1, etc.

    For the 2D index layout, see EdgeGrid class docstring. The corresponding
    1D index layout is::

        '        /           /           /           /           /
        '      42          44          46           48          50
        '      /           /           /           /           /
        '     +-----43----+-----45----+-----47----+-----49----+
        '    /|          /|          /|          /|          /|
        '   / |         / |         / |         / |         / |
        '    37  /       38  /       39  /       40  /       41  /
        '     | 28        | 30        | 32        | 34        | 36
        '     |/          |/          |/          |/          |/
        '     +-----29----+-----31----+-----33----+-----35----+
        '    /|          /|          /|          /|          /|
        '   / |         / |         / |         / |         / |
        '    23  /       24  /       25  /       26  /       27  /
        '     | 14        | 16        | 18        | 20        | 22
        '     |/          |/          |/          |/          |/
        '     +-----15----+-----17----+-----19----+-----21----+
        '    /|          /|          /|          /|          /|
        '   / |         / |         / |         / |         / |
        '     9  /       10  /       11  /       12  /       13  /
        '     | 0         | 2         | 4         | 6         | 8
        '     |/          |/          |/          |/          |/
        '     +-----1-----+-----3-----+-----5-----+-----7-----+
        '    /           /           /           /           /
        '   /           /           /           /           /

    :param index: The 0-indexed edge index in the Fortran 1D system
    :type index: int
    :param rows: The number of rows in the grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the grid
    :type nodes_in_row: int
    :return: The (row, rank) coordinates for the edge in the 2D system
    :rtype: Tuple[int, int]
    :raises IndexError: If index is out of bounds [0, total_edges-1]
    """
    # The number of edges in a full row: 3 of each per node, except the last
    # node which has no x-edge.
    full_row = 3 * nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except
    # the last node which has no x-edge.
    xz_row = 2 * nodes_in_row - 1
    # The total number of edges in the grid: full_row for each row, except the
    # last row which has no y-edges.
    total_edges = full_row * (rows - 1) + xz_row

    # Check if the index given is in-bounds.
    if index < 0 or index > total_edges - 1:
        raise IndexError(
            f"Index ({index}) out of bounds. Edge indices are in the range "
            f"[0..{total_edges-1}]."
        )

    # Count the number of full rows before this edge
    i = index // full_row
    # Determine the number of edges (in 1-D order) before this one in its own row
    rank = index % full_row
    # If all x- and z-edges are already counted, this must be a y-edge
    if rank > xz_row - 1:
        # Its index in the list of y-edges is the index of its triplet in the 2-D index
        triplet = rank - xz_row
        # The y-fiber is first in its triplet, so count up the triplets before this one.
        j = triplet * 3
    else:
        # Else it is an x- or z-edge. So find out which triplet it is in by
        # counting pairs of x- and z-edges.
        triplet = rank // 2
        # Then we need to insert all the y-edges for the preceding triplets,
        # and the y-edge for this triplet.
        j = rank + triplet + 1

    # Return the co-ordinates in the 2-D ordering.
    return i, j


def from_fortran_edge_index_array(
    index_array: np.ndarray, rows: int, nodes_in_row: int
) -> np.ndarray:
    """Convert array of Fortran 1D edge indices to 2D (i, j) coordinates (vectorized).

    Vectorized version of from_fortran_edge_index() that operates on entire arrays
    at once, providing significant performance improvements for bulk conversions.

    Converts multiple Fortran 1D indices to their corresponding (row, rank) pairs
    in the Python 2D indexing system.

    :param index_array: 1D array of 0-indexed edge indices in the Fortran system.
        Shape: (n_edges,)
    :type index_array: np.ndarray
    :param rows: The number of rows in the grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the grid
    :type nodes_in_row: int
    :return: 2D array where result[k, :] = (i, j) coordinates for index_array[k].
        Shape: (n_edges, 2), dtype: uint32
    :rtype: np.ndarray
    :raises IndexError: If any indices are out of bounds [0, total_edges-1]
    """
    # The number of edges in a full row: 3 of each per node, except the last
    # node which has no x-edge.
    full_row = 3 * nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except
    # the last node which has no x-edge.
    xz_row = 2 * nodes_in_row - 1
    # The total number of edges in the grid: full_row for each row, except the
    # last row which has no y-edges.
    total_edges = full_row * (rows - 1) + xz_row

    # Check if the index given is in-bounds.
    if np.count_nonzero(index_array < 0) or np.count_nonzero(
        index_array > total_edges - 1
    ):
        raise IndexError(
            f"Index ({index_array}) out of bounds. Edge indices are in the range "
            f"[0..{total_edges-1}]."
        )

    out_array = np.empty(shape=(index_array.size, 2), dtype=np.uint32)
    # Count the number of full rows before this edge
    # Determine the number of edges (in 1-D order) before this one in its own row
    out_array[:, 0], rank = np.divmod(index_array, full_row)
    # If all x- and z-edges are already counted, this must be a y-edge
    # Its index in the list of y-edges is the index of its triplet in the 2-D index
    # The y-fiber is first in its triplet, so count up the triplets before this one.
    out_array[rank > xz_row - 1, 1] = (rank[rank > xz_row - 1] - xz_row) * 3
    # Else it is an x- or z-edge. So find out which triplet it is in by
    # counting pairs of x- and z-edges.
    # Then we need to insert all the y-edges for the preceding triplets,
    # and the y-edge for this triplet.
    out_array[rank <= xz_row - 1, 1] = (
        rank[rank <= xz_row - 1] + rank[rank <= xz_row - 1] // 2 + 1
    )
    # Return the co-ordinates in the 2-D ordering.
    return out_array


def to_fortran_edge_index(i: int, j: int, rows: int, nodes_in_row: int) -> int:
    """Convert 2-dimensional (i, j) coordinates to a 1-dimensional Fortran edge index.

    Converts from the (i, j) index used by this Python package to the 1D index
    used in Fortran "Macro" data structures and files.

    .. note::
        **All indices are zero-indexed.** Add 1 before using in Fortran code.

    This function aids compatibility with legacy Fortran code during the transition
    phase and should be deprecated once Python implementation is complete.

    For the 2D index layout, see EdgeGrid class docstring. For the 1D layout,
    see from_fortran_edge_index() docstring.

    :param i: The row index of the edge (0 to rows-1)
    :type i: int
    :param j: The rank index of the edge within its row (0 to full_row-1)
    :type j: int
    :param rows: The number of rows in the grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the grid
    :type nodes_in_row: int
    :return: The 0-indexed edge address in the Fortran Macro structure
    :rtype: int
    :raises IndexError: If i or j are out of bounds, or if (i,j) refers to a
        y-edge on the top row (which doesn't exist)
    """
    # The number of edges in a full row: 3 of each per node, except the last
    # node which has no x-edge.
    full_row = 3 * nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except
    # the last node which has no x-edge.
    xz_row = 2 * nodes_in_row - 1

    # Check that the indices given are valid
    if i < 0 or i > rows - 1:
        raise IndexError(f"Index i={i} out of bounds. Rows are [0..{rows-1}].")
    if j < 0 or j > full_row - 1:
        raise IndexError(
            f"Index j={j} out of bounds. Edges in each row are [0..{full_row-1}]."
        )

    # Add up the number of edges that are in the preceding rows.
    index = i * full_row
    # Determine which y-, z-, and x-edge triplet in its row it belongs to.
    triplet = j // 3
    if j % 3 == 1:
        # If the edge is a z-edge,
        # the number of edges before it in this row is two (x- and z-edges) per triplet.
        index += 2 * triplet
    elif j % 3 == 2:
        # If the edge is an x-edge,
        # the number of edges before it in this row is two (x- and z-edges) per triplet,
        # plus it's partner z-edge.
        index += 2 * triplet + 1
    elif j % 3 == 0:
        # If the edge is a y-edge
        if i < rows - 1:
            # In all but the last row,
            # the number of edges before it is all x- and z-edges in its row,
            # plus the y-edges before it.
            index += xz_row + triplet
        else:
            # In the last row, there are no y-edges.
            raise IndexError(f"No y-edges on the top row (row {rows-1}).")

    return index


def to_fortran_edge_index_array(
    index_array: np.ndarray, rows: int, nodes_in_row: int
) -> np.ndarray:
    """Convert array of 2D (i, j) coordinates to Fortran 1D edge indices (vectorized).

    Vectorized version of to_fortran_edge_index() that operates on entire arrays
    at once, providing significant performance improvements for bulk conversions.

    Converts multiple (row, rank) coordinate pairs from the Python 2D system to
    their corresponding Fortran 1D indices.

    :param index_array: 2D array where each row is (i, j) coordinates in Python system.
        Shape: (n_edges, 2)
    :type index_array: np.ndarray
    :param rows: The number of rows in the grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the grid
    :type nodes_in_row: int
    :return: 1D array of 0-indexed Fortran edge indices corresponding to input coordinates.
        Shape: (n_edges,), dtype: uint32. Add 1 for Fortran code.
    :rtype: np.ndarray
    :raises IndexError: If any coordinates are out of bounds or refer to non-existent
        y-edges on the top row
    """
    # The number of edges in a full row: 3 of each per node, except the last
    # node which has no x-edge.
    full_row = 3 * nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except
    # the last node which has no x-edge.
    xz_row = 2 * nodes_in_row - 1

    # Check that the indices given are valid
    if np.count_nonzero(index_array[:, 0] < 0) or np.count_nonzero(
        index_array[:, 0] > rows - 1
    ):
        raise IndexError(
            f"Index i={index_array[:, 0]} out of bounds. Rows are [0..{rows-1}]."
        )
    if np.count_nonzero(index_array[:, 1] < 0) or np.count_nonzero(
        index_array[:, 1] > full_row - 1
    ):
        raise IndexError(
            f"Index j={index_array[:, 1]} out of bounds. Edges in each row are [0..{full_row-1}]."
        )
    # Determine which y-, z-, and x-edge triplet in its row it belongs to.
    triplet, direction = np.divmod(index_array[:, 1], 3)
    # In the last row, there are no y-edges.
    if np.count_nonzero(index_array[direction == 0, 0] >= rows - 1):
        raise IndexError(f"No y-edges on the top row (row {rows-1}).")
    out = np.empty(shape=index_array.shape[0], dtype=np.uint32)
    # Add up the number of edges that are in the preceding rows.
    out = index_array[:, 0] * full_row
    # If the edge is a z-edge,
    # the number of edges before it in this row is two (x- and z-edges) per triplet.
    out[direction == 1] += 2 * triplet[direction == 1]
    # If the edge is an x-edge,
    # the number of edges before it in this row is two (x- and z-edges) per triplet,
    # plus its partner z-edge.
    out[direction == 2] += 2 * triplet[direction == 2] + 1
    # If the edge is a y-edge
    # In all but the last row,
    # the number of edges before it is all x- and z-edges in its row,
    # plus the y-edges before it.
    out[direction == 0] += xz_row + triplet[direction == 0]

    return out
