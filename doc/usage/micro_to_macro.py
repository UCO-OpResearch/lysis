#!/usr/bin/env python3

import argparse
import os

import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_code", type=str)
    parser.add_argument(
        "--in_code",
        type=str,
        default="",
        help="The code to add to all input filenames (should include any leading underscores, but NOT the file extension).",
    )
    parser.add_argument(
        "--out_code",
        type=str,
        default="",
        help="The code to add to all output filenames (should include any leading underscores, but NOT the file extension).",
    )
    parser.add_argument(
        "-F", "--rows",
        type=int,
        default="",
        help="The number of rows in the macroscale simulation EdgeGrid. This should include any non-fiber rows."
    )
    parser.add_argument(
        "-N", "--cols",
        type=int,
        default="",
        help="The number of nodes in one row of the macroscale simulation EdgeGrid."
    )
    return parser.parse_args()

class Neighbors:
    def __init__(self):
        self.X = ((-1, -1, 0, 0, 0, 0, 0, 0), (-2, 1, -2, 1, -1, -1, 2, 2))
        self.Y = ((0, 0, 1, 1, 0, 1, 0, 1), (1, 1, 1, 1, -1, -1, 2, 2))
        self.Z = ((-1, -1, 0, 0, 0, 0, 0, 0), (-1, -1, -1, -1, -2, -2, 1, 1))
        self.TOP_REFL = (0, 0, -1, -1, 0, 0, 0, 0)
        self.BOTTOM_REFL = (1, 1, 0, 0, 0, 0, 0, 0)
        self.LEFT_REFL = (0, 0, 0, 0, 3, 3, 0, 0)
        self.RIGHT_REFL = (0, 0, 0, 0, 0, 0, -3, -3)


def full_row(rows: int, nodes_in_row: int) -> int:
    """
    Calculates the number of edges in a full row of the edge grid

    :param rows: The number of rows in the edge grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the edge grid
    :type nodes_in_row: int
    :return: The number of edges in a full row of the edge grid
    :rtype: int
    """
    return 3 * nodes_in_row - 1

def xz_row(rows: int, nodes_in_row: int) -> int:
    """
    Calculates the number of x- and z-edges in a row of the edge grid

    :param rows: The number of rows in the edge grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the edge grid
    :type nodes_in_row: int
    :return: The number of x- and z-edges in a row of the edge grid
    :rtype: int
    """
    return 2 * nodes_in_row - 1

def total_edges(rows: int, nodes_in_row: int) -> int:
    """
    Calculates the total number of edges in an edge grid

    :param rows: The number of rows in the edge grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the edge grid
    :type nodes_in_row: int
    :return: The total number of edges in an edge grid
    :rtype: int
    """
    return full_row(rows, nodes_in_row) * (rows - 1) + xz_row(
        rows, nodes_in_row
    )

def from_fortran_edge_index_array(
    index_array: np.ndarray, rows: int, nodes_in_row: int
) -> np.ndarray:
    """
    Does the same as the from_fortran_edge_index() method, but with a whole numpy array at once.

    :param index_array: A (-1,) array of 0-indexed indices in the Fortran 1-D indexing system.
    :type index_array: np.ndarray
    :param rows: The number of rows in the grid
    :type rows: int
    :param nodes_in_row: The number of nodes each row of the grid
    :type nodes_in_row: int
    :raises IndexError: Raised if any of the indices are out of bounds
    :return: A (-1, 2) array with each row containing the indices of an edge in the 2-D ordering system
    :rtype: np.ndarray
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


def to_fortran_edge_index_array(
    index_array: np.ndarray, rows: int, nodes_in_row: int
) -> np.ndarray:
    """
    Does the same as the to_fortran_edge_index() method, but with a whole numpy array at once.

    :param index_array: A (-1, 2) array with each row containing the indices of an edge in the 2-D ordering system
    :type index_array: np.ndarray
    :param rows: The number of rows in the grid.
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the grid.
    :type nodes_in_row: int
    :raises IndexError: Raised if any of the indices are out of bounds
    :return: A (-1,) array of 0-indexed indices in the Fortran 1-D indexing system.
    :rtype: np.ndarray
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


def generate_fortran_neighborhood_structure(rows: int, nodes_in_row: int) -> np.ndarray:
    """
    Generates the neighbor structure needed by the Fortran Macroscale code.
    This is a (-1, 8) array with a row for each edge, ordered by their Fortran index.
    Each row contains the 0-indexed, fortran (1-D) indices for the 8 neighbors of the edge.
    Each row is sorted in increasing order for consistency with the Fortran code.

    :param rows: The number of rows in the edge grid
    :type rows: int
    :param nodes_in_row: The number of nodes in each row of the edge grid
    :type nodes_in_row: int
    :return: A (-1, 8) NumPy array of dtype uint32
    :rtype: np.ndarray
    """
    NEIGHBORHOOD = Neighbors()
    # Generate a list of all fortran-index edges and get the equivalent 2-d indeces
    edges = from_fortran_edge_index_array(
        np.arange(total_edges(rows, nodes_in_row)), rows, nodes_in_row
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
    right = edges[:, 1] >= full_row(rows, nodes_in_row) - 2
    # Add an axis in the middle
    edges = edges.reshape(-1, 1, 2)
    # and then duplicate 8 times along that axis.
    # This is where the 8 neighbors will go
    neighbors = np.repeat(edges, 8, axis=1)

    # Create a matrix for the neighbor delta
    adds = np.empty(neighbors.shape, dtype=int)
    # Get the appropriate neighbor delta for each edge type (x, y, z) and put them in the matrix
    adds[edge_type == 0] = np.array(NEIGHBORHOOD.Y).T
    adds[edge_type == 1] = np.array(NEIGHBORHOOD.Z).T
    adds[edge_type == 2] = np.array(NEIGHBORHOOD.X).T
    # Deal with the boundaries
    adds[top, :, 0] += NEIGHBORHOOD.TOP_REFL
    adds[bottom, :, 0] += NEIGHBORHOOD.BOTTOM_REFL
    adds[left, :, 1] += NEIGHBORHOOD.LEFT_REFL
    adds[right, :, 1] += NEIGHBORHOOD.RIGHT_REFL
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


# Get command-line arguments
args = parse_arguments()

# Create neighborhood structure and save to file
fort_neighbors = generate_fortran_neighborhood_structure(args.rows, args.cols) + 1
fort_neighbors.tofile(os.path.join("data", args.run_code, "neighbors.dat"), sep=os.linesep)

# Read in microscale data
tpa_leaving_time = np.fromfile(os.path.join("data", args.in_code, f"tPA_time{args.in_code}.dat"))
fiber_degraded = np.fromfile(
    os.path.join("data", args.in_code, f"lyscomplete{args.in_code}.dat"), 
    dtype=np.int32,
).astype(bool)
sim_final_time = np.fromfile(os.path.join("data", args.in_code, f"lysis{args.in_code}.dat"))



# Get the number of microscale runs and set the dimensions of the bins so that we get 100 bins
set_size = tpa_leaving_time.size // 100
# tPAleave is the CDF of the tPA leaving time distribution.
# This is really just a list of edgepoints from the bins for tPA leaving time
# These bins are evenly distributed along the interval [0, 1]
tPAleave = np.append(np.arange(0, 1, 0.01), [1.0])
np.savetxt(os.path.join("data", args.run_code, f"tPAleave{args.in_code}.dat"), tPAleave)

# The remaining data will be arranged into 100 bins
# according to the time tPA left the simulation.
# Get the sorted ordering of the tPA leaving times
indices = tpa_leaving_time.argsort()
# Find the tPA leaving times for the edges of each bin.
tsectPA = np.append([0], tpa_leaving_time[indices[set_size - 1 :: set_size]])
np.savetxt(os.path.join("data", args.run_code, f"tsectPA{args.in_code}.dat"), tsectPA)

# If full degradation did NOT occur,
# this matrix currently contains the ending time of the simulation.
# Replace these times with an 'infinity' marker of 6,000 seconds
sim_final_time[~fiber_degraded] = 6000
# Rearrange the matrix so that each row contains the lysis times for simulations
# corresponding to the matching bin in the ``bin_edge_proportions`` vector.
# Then sort the rows (bins) individually by lysis time.
# Finally, transpose the matrix so that the bins are arranged in columns.
lysismat = np.stack(
    [np.sort(sim_final_time[indices[i * set_size : (i + 1) * set_size]]) for i in range(100)]
).T
np.savetxt(os.path.join("data", args.run_code, f"lysismat{args.in_code}.dat"), lysismat)

# Find the location of the first '6000' entry in each column of the ``binned_fiber_degrade_time`` matrix
# Then convert to 1-indexing.
lenlysisvect = lysismat.argmax(axis=0)+1
np.savetxt(os.path.join("data", args.run_code, f"lenlysisvect{args.in_code}.dat"), lenlysisvect)
