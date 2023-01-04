from typing import Any, Callable, List, Tuple

import numpy as np

from .util import Const, BoundaryCondition, Experiment


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

CONST = Const()


class EdgeGrid(object):
    """The main class containing a 3-D grid of edges. This represents an xy-planar slice, one edge high, of a clot.

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
    def __init__(self,
                 exp: Experiment,
                 thisgrid_rows: int | None = None,
                 thisgrid_starting_row: int = 0,
                 boundary_conditions: Tuple[BoundaryCondition, BoundaryCondition] | None = None,
                 initial_fiber_status: float = float('inf'),
                 ):
        """Initializes an EdgeGrid.

        Args:
            exp: The experiment that this EdgeGrid is a part of.
                This structure passes many of the parameters used to set up the experiment.
            thisgrid_rows: The number of rows in this EdgeGrid
                **NOTE**: This does not count duplicated edges used with CONTINUING boundary conditions.
                That is, if this EdgeGrid is in the middle of a larger grid,
                there will actually be two more rows in this grid,
                "Row thisgrid_starting_row-1", which will represent the top row of the grid below, and
                "Row thisgrid_starting_row + thisgrid_rows", which will represent the bottom row of the grid above.
            thisgrid_starting_row: The first row in this EdgeGrid (in the larger grid numbering system.
            boundary_conditions: A tuple that maps boundaries from CONST.BOUND to conditions from CONST.BOUND_COND.
            initial_fiber_status: The initial value of the fiber status.
                This should generally be an infinite degrade time (degrade time > experiment length).
        """
        # TODO(bpaynter): Redo so that externally the row numbers represent the numbering in the larger grid
        #                   while still starting at zero internally.
        #                   That is, if there are 100 rows total, split into 5 grids (0 through 4),
        #                   of which I am Grid[3], then externally my rows are 60 through 79.
        #                   Internally I have storage for 22 rows which represent global rows 59 through 80.
        #   2023-01-03:  Was busy with this but need to finish!

        self.total_rows = exp.macro_params.rows
        """int: The total number of rows in the entire grid (of which this EdgeGrid is either a part or the whole)."""
        self.nodes_in_row = exp.macro_params.cols
        """int: The number of nodes (not edges) in each row of this EdgeGrid."""
        self.total_empty_rows = exp.macro_params.empty_rows
        """int: The number of empty (fibrin-free) rows in the entire grid.
                There are only empty rows in this EdgeGrid if total_empty_rows > thisgrid_starting_row.
                These will always be the bottom rows of the grid (i.e., "Row 0" through "Row empty_rows-1")"""
        self.thisgrid_rows = thisgrid_rows if thisgrid_rows is not None else self.total_rows
        self.thisgrid_starting_row = thisgrid_starting_row
        self.thisgrid_empty_rows = max(self.total_empty_rows - self.thisgrid_starting_row, 0)
        """int: Number of empty rows in this EdgeGrid. If they exist they are 
                "Row thisgrid_starting_row" through "Row max(thisgrid_starting_row, total_empty_rows)-1"
                If a top and/or bottom "shadow" row exists, it is NOT counted in this number. 
                This row is always treated as empty because, even if it does have fibrin, 
                those fibers exist in a different EdgeGrid."""
        # Sanity check
        if self.thisgrid_starting_row + self.thisgrid_rows > self.total_rows:
            raise AttributeError('This EdgeGrid overruns the top of the total grid.')
        # Set the appropriate boundary conditions
        if boundary_conditions is not None:
            self.boundary_conditions = boundary_conditions
        elif self.thisgrid_rows == self.total_rows:                                 # We are the entire grid
            self.boundary_conditions = (CONST.BOUND_COND.REFLECTING, CONST.BOUND_COND.REFLECTING)
        elif self.thisgrid_starting_row == 0:                                       # We are the bottom EdgeGrid
            self.boundary_conditions = (CONST.BOUND_COND.CONTINUING, CONST.BOUND_COND.REFLECTING)
        elif self.thisgrid_rows + self.thisgrid_starting_row == self.total_rows:    # We are the top EdgeGrid
            self.boundary_conditions = (CONST.BOUND_COND.REFLECTING, CONST.BOUND_COND.CONTINUING)
        else:                                                                       # We are a middle EdgeGrid
            self.boundary_conditions = (CONST.BOUND_COND.CONTINUING, CONST.BOUND_COND.CONTINUING)

        self.edges_in_row = 3*self.nodes_in_row - 1
        self.fiber_rows = self.thisgrid_rows - self.thisgrid_empty_rows
        self._fiber_status = initial_fiber_status * np.ones((self.fiber_rows, self.edges_in_row),
                                                            dtype=np.double)

        row_shift: Callable[[int], int] = lambda i: i - self.thisgrid_starting_row - self.thisgrid_empty_rows
        identity: Callable[[int], int] = lambda j: j
        self.fiber_status = self._ShiftedArrayAccess(self._fiber_status, (row_shift, identity))
        """np.ndarray: The status of the fibers in this EdgeGrid. This is essentially the degrade time of the fiber.
            If degrade time < current time, the fiber is degraded.
            e.g., self.fiber_status[60, 25] is the status of the fiber located along edge [60, 25].
            
            Supports two-dimensional slicing.
            
            Note that the location should be given in reference to the entire grid, not just this EdgeGrid."""

        self.molecules = None
        """List[List[List]]: A structure containing the indices of the molecules currently
            located at each edge of this EdgeGrid. 
            e.g., self.molecules[60, 25] is a list of molecules currently located at edge [60, 25].
             
            Supports slicing in one dimension only.
             
            Note that the location should be given in reference to the entire grid, not just this EdgeGrid."""
        self._molecules = [[[] for _ in range(self.edges_in_row)] for _ in range(self.thisgrid_rows)]
        if self.boundary_conditions[CONST.BOUND.TOP] == CONST.BOUND_COND.CONTINUING:
            self._molecules.append([[] for _ in range(self.edges_in_row)])
        if self.boundary_conditions[CONST.BOUND.BOTTOM] == CONST.BOUND_COND.CONTINUING:
            self._molecules.append([[] for _ in range(self.edges_in_row)])

            row_shift: Callable[[int], int] = lambda i: i - self.thisgrid_starting_row + 1
            identity: Callable[[int], int] = lambda j: j
            self.molecules = self._ShiftedArrayAccess(self._molecules, (row_shift, identity))
        else:
            row_shift: Callable[[int], int] = lambda i: i - self.thisgrid_starting_row
            identity: Callable[[int], int] = lambda j: j
            self.molecules = self._ShiftedArrayAccess(self._molecules, (row_shift, identity))

    def _is_valid_index(self, i: int, j: int) -> str | None:
        """Checks if [i, j] is a valid index for this EdgeGrid.

        Args:
            i: The index of the edge's row.
            j: The index of the edge within its row.

        Returns: None if the index is valid, or an error message (as a string) if the index is invalid.
        """
        # Edges are 0 through self.edges_in_row-1
        if j < 0 or j > self.edges_in_row - 1:
            return (f'Index j={j} out of bounds. '
                    f'This model only has edges [0..{self.edges_in_row - 1}] in each row.')
        # Rows are 0 through self.rows-1
        if i < 0 or i > self.thisgrid_rows - 1:
            return f'Index i={i} out of bounds. This model only has rows [0..{self.thisgrid_rows - 1}].'
        # If this EdgeGrid is at the top of a larger sliced grid, or is the whole grid itself,
        # then the top row has no y-edges in it.
        if (
                i == self.thisgrid_rows - 1
                and
                self.boundary_conditions[CONST.BOUND.TOP] == CONST.BOUND_COND.REFLECTING
                and
                j % 3 == 0
           ):
            return f'y-edges do not exist on the top row of this grid (row {self.thisgrid_rows - 1}). Location ({i}, {j})'
        # Everything seems fine, so return None
        return None

    def neighbor(self, i: int, j: int, k: int) -> Tuple[int, int]:
        """Finds the co-ordinates of a neighboring edge.

        Given the co-ordinates of an edge, and the index (in the neighborhood) of a neighboring edge,
        this method will find the co-ordinates of the requested neighbor,
        taking appropriate boundary conditions into account.

        Neighborhoods each contain eight other edges and are constructed as in the following examples
        (without boundary conditions). Note that in each situation, the generating edge is labelled "(i,j)"
        and the members of the neighborhood are labelled "k:(i+p,j+q)" where (p,q) is the vector from the
        generating edge to the neighbor.
        ::

            * x-edge neighborhood::
                
                '        2:(i,j-2)     3:(i,j+1)    
                '            |             |    
                '            |  /          |  /         
                '            | 5:(i,j-1)   | 7:(i,j+2)  
                '            |/            |/         
                '            +----(i,j)----+           
                '          / |            /|           
                ' 4:(i,j-1)  |   6:(i,j+2) |           
                '        /   |          /  |     
                '            |             |      
                '      0:(i-1,j-2)   1:(i-1,j+1)       

            * y-edge neighborhood::

                '                    /
                '                   3:(i+1,j+1)    
                '                  /             
                ' 5:(i+1,j-1)-----+-----7:(i+1,j+2)  
                '                /|               
                '     2:(i+1,j+1) |               
                '              /  |               
                '               (i,j)             
                '                 |  /           
                '                 | 1:(i,j+1)      
                '                 |/             
                '   4:(i,j-1)-----+-----6:(i,j+2)    
                '                /                
                '       0:(i,j+1)  
                '              /           

            * z-edge neighborhood::      

                '                         |           
                '               3:(i+1,j-1)       
                '                         |         
                '         5:(i+1,j-2)-----+-----7:(i+1,j+1)    
                '                        /|            
                '                       / 1:(i,j-1)          
                '                      /  |
                '                   (i,j) 
                '                 |  /                  
                '       2:(i+1,j-1) /                   
                '                 |/                  
                ' 4:(i+1,j-2)-----+-----6:(i+1,j+1)                 
                '                 |        
                '                 0:(i,j-1)        
                '                 |  

        Args:
            i: The row of the generating edge.
            j: The index of the generating edge within its row.
            k: The index of the neighbor in the neighborhood (i.e., the (k+1)st neighbor)

        Returns: The pair of co-ordinates (in the EdgeGrid) of the requested neighbor.
        """
        # Check that we are in-bounds

        # Neighbors are 0 through 7
        if k < 0 or k > 7:
            raise IndexError(f'Index k={k} out of bounds. Each neighborhood has items [0..7].')
        # Check if i and j are valid
        valid_index = self._is_valid_index(i, j)
        if valid_index is not None:
            raise IndexError(valid_index)

        # The index of the neighboring fiber being requested
        neighbor_i = i
        neighbor_j = j

        # Move to the required neighbor
        if j % 3 == 0:              # We are a y-edge
            neighbor_i += CONST.NEIGHBORHOOD.Y[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.Y[1][k]
        elif j % 3 == 1:            # We are a z-edge
            neighbor_i += CONST.NEIGHBORHOOD.Z[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.Z[1][k]
        elif j % 3 == 2:            # We are an x-edge
            neighbor_i += CONST.NEIGHBORHOOD.X[0][k]
            neighbor_j += CONST.NEIGHBORHOOD.X[1][k]

        # Deal with boundary conditions

        # The bottom boundary of the grid.
        # Note that, if the edge generating the neighborhood is a y-edge,
        # then its neighborhood only involves fibers on its own row, or the row above.
        # Thus, the neighborhood of y-edges never overruns the bottom of the grid
        if (
                i == 0                          # We are at the bottom of the grid
                and
                j % 3 > 0                       # and it is a z- or x-edge
           ):
            if self.boundary_conditions[CONST.BOUND.BOTTOM] == CONST.BOUND_COND.REFLECTING:
                # If the bottom is reflective, then we simply shift any neighbors on row -1, to the row above.
                neighbor_i += CONST.NEIGHBORHOOD.BOTTOM_REFL[k]

        # The top boundary of the grid.
        # Note that, if the edge generating the neighborhood is a y-edge,
        # then its neighborhood only involves fibers on its own row, or the row above.
        # But, if this row is the top of the entire experiment (REFLECTING) then there are no y-edges on this row
        # If this row is the top of one slice (CONTINUING) then this row represents the bottom row of the next slice
        # and should not be processed here
        elif (
                i == self.thisgrid_rows-1                # We are at the top of the grid
                and
                j % 3 > 0                       # and it is a z- or x-edge
             ):
            if self.boundary_conditions[CONST.BOUND.TOP] == CONST.BOUND_COND.REFLECTING:
                neighbor_i += CONST.NEIGHBORHOOD.TOP_REFL[k]

        # The left boundary of the grid.
        # Note that, if the edge generating the neighborhood is an x-edge,
        # Then its neighborhood never overruns the side of the grid
        if j <= 1:                              # We are the left-most y- or z-edge
            neighbor_j += CONST.NEIGHBORHOOD.LEFT_REFL[k]

        # The right boundary of the grid.
        # Note that, if the edge generating the neighborhood is an x-edge,
        # Then its neighborhood never overruns the side of the grid
        elif j >= self.edges_in_row-2:          # We are the right-most y- or z-edge
            neighbor_j += CONST.NEIGHBORHOOD.RIGHT_REFL[k]

        # Return the co-ordinates of the requested neighbor.
        return neighbor_i, neighbor_j

    def pop_molecule_row(self, i: int) -> List[List[int]]:
        if -1 < i < self.thisgrid_rows:
            raise IndexError('You should not be removing molecules from the actual rows of this grid. '
                             'Only shadow rows allowed.')
        row = self.molecules[i].copy()
        self.molecules[i] = [[] for _ in range(self.edges_in_row)]
        return row

    def place_molecule_row(self, i: int, row: List[List[int]]):
        if 0 < i < self.thisgrid_rows-1:
            raise IndexError('You should only add a row to the first or last rows of this grid.')
        for j in range(self.edges_in_row):
            self.molecules[i, j].append(row[j])

    def move_molecule(self, m: int, start: Tuple[int, int], end: Tuple[int, int]):
        self.molecules[start].remove(m)
        self.molecules[end].append(m)


    class _ShiftedArrayAccess:
        def __init__(self,
                     array: List[List[Any]],
                     shift: Tuple[Callable[[int], int], ...],
                     ):
            self.array = array
            self.shift = shift

        @staticmethod
        def _shift_key(key: slice | int,
                       shift_func: Callable[[int], int]
                       ) -> slice | int:
            if isinstance(key, slice):
                start = shift_func(key.start) if key.start is not None else None
                stop = shift_func(key.stop) if key.stop is not None else None
                return slice(start, stop, key.step)
            else:
                return shift_func(key)

        def __getitem__(self, key: int | slice | Tuple[int | slice, ...]) -> Any:
            if isinstance(key, Tuple):
                shifted_keys = tuple([self._shift_key(k, self.shift[idx]) for idx, k in enumerate(key)])
            else:
                shifted_keys = (self._shift_key(key, self.shift[0]),)
            if isinstance(self.array, np.ndarray):
                return self.array[shifted_keys]
            else:
                out = self.array
                for k in shifted_keys:
                    out = out[k]
                return out

        def __setitem__(self, key: int | slice | Tuple[int | slice], value: Any):
            if isinstance(key, Tuple):
                shifted_keys = tuple([self._shift_key(k, self.shift[idx]) for idx, k in enumerate(key)])
            else:
                shifted_keys = (self._shift_key(key, self.shift[0]),)
            if isinstance(self.array, np.ndarray):
                self.array[shifted_keys] = value
            else:
                place = self.array
                for k in shifted_keys[:-1]:
                    place = place[k]
                place[shifted_keys[-1]] = value


def from_fortran_edge_index(index: int, rows: int, nodes_in_row: int) -> Tuple[int, int]:
    """Converts a 1-dimensional index of an edge to its 2-dimensional index.

    Converts from the index for an edge used in the "Macro" Fortran data structures and data files
    into the (i, j) index used by this package.

    **NOTE**: All indices in this method are zero-indexed.
    They will need to be incremented by one if used in Fortran code.

    This method exists to aid with compatability with the Fortran "Macro" code during the transition phase.
    It should be depreciated once the Python code is complete.

    For an example of the 2-dimensional index, see the DocString for the EdgeGrid class.
    The 1-dimensional index arrangement for the same setup is as follows::

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

    Args:
        index: The index of an edge in the 1-dimensional system
        rows: The number of rows in the grid
        nodes_in_row: The number of nodes each row of the grid

    Returns: The pair of co-ordinates for the edge in the 2-dimensional system.
    """
    # The number of edges in a full row: 3 of each per node, except the last node which has no x-edge.
    full_row = 3*nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except the last node which has no x-edge.
    xz_row = 2*nodes_in_row - 1
    # The total number of edges in the grid: full_row for each row, except the last row which has no y-edges.
    total_edges = full_row * (rows-1) + xz_row

    # Check if the index given is in-bounds.
    if index < 0 or index > total_edges-1:
        raise IndexError(f'Index ({index}) out of bounds. Edge indices are in the range [0..{total_edges-1}].')

    # Count the number of full rows before this edge
    i = index // full_row
    # Determine the number of edges (in 1-D order) before this one in its own row
    index_in_row = index % full_row
    # If all x- and z-edges are already counted, this must be a y-edge
    if index_in_row > xz_row-1:
        # Its index in the list of y-edges is the index of its triplet in the 2-D index
        triplet = index_in_row - xz_row
        # The y-fiber is first in its triplet, so count up the triplets before this one.
        j = triplet * 3
    else:
        # Else it is an x- or z-edge. So find out which triplet it is in by counting pairs of x- and z-edges.
        triplet = index_in_row // 2
        # Then we need to insert all the y-edges for the preceding triplets, and the y-edge for this triplet.
        j = index_in_row + triplet + 1

    # Return the co-ordinates in the 2-D ordering.
    return i, j


def to_fortran_edge_index(i: int, j: int, rows: int, nodes_in_row: int) -> int:
    """Converts a 2-dimensional index of an edge to its 1-dimensional index.

    Converts from the (i, j) index for an edge used by this package into the index used in
    "Macro" Fortran data structures and data files.

    **NOTE**: All indices in this method are zero-indexed.
    They will need to be incremented by one if used in Fortran code.

    This method exists to aid with compatability with the Fortran "Macro" code during the transition phase.
    It should be depreciated once the Python code is complete.

    For an example of the 2-dimensional index, see the DocString for the EdgeGrid class.
    For a matching example of the 1-dimensional index, see the DocString for the from_fortran_edge_index() method.

    Args:
        i: The index of the edge's row.
        j: The index of the edge within its row.
        rows: The number of rows in the grid.
        nodes_in_row: The number of nodes in each row of the grid.

    Returns: A zero-indexed address of the edge in the Fortran Macro structure.
    """
    # The number of edges in a full row: 3 of each per node, except the last node which has no x-edge.
    full_row = 3*nodes_in_row - 1
    # The number of x- and z-edges in a full row: 2 of each per node, except the last node which has no x-edge.
    xz_row = 2*nodes_in_row - 1

    # Check that the indices given are valid
    if i < 0 or i > rows-1:
        raise IndexError(f'Index i={i} out of bounds. Rows are [0..{rows-1}].')
    if j < 0 or j > full_row-1:
        raise IndexError(f'Index j={j} out of bounds. Edges in each row are [0..{full_row-1}].')

    # Add up the number of edges that are in the preceding rows.
    index = i * full_row
    # Determine which y-, z-, and x-edge triplet in its row it belongs to.
    triplet = j // 3
    if j % 3 == 1:
        # If the edge is a z-edge,
        # the number of edges before it in this row is two (x- and z-edges) per triplet.
        index += 2*triplet
    elif j % 3 == 2:
        # If the edge is an x-edge,
        # the number of edges before it in this row is two (x- and z-edges) per triplet,
        # plus it's partner z-edge.
        index += 2*triplet + 1
    elif j % 3 == 0:
        # If the edge is a y-edge
        if i < rows-1:
            # In all but the last row,
            # the number of edges before it is all x- and z-edges in its row,
            # plus the y-edges before it.
            index += xz_row + triplet
        else:
            # In the last row, there are no y-edges.
            raise IndexError(f'No y-edges on the top row (row {rows-1}).')

    return index
