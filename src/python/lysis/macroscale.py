from typing import Tuple

import numpy as np
import util

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

CONST = util.Const()


class FiberGrid(object):
    def __init__(self,
                 nodes_in_col: int,
                 nodes_in_row: int,
                 empty_rows: int,
                 boundary_conditions: Tuple[util.BoundaryCondition,
                                            util.BoundaryCondition]
                 ):
        self.nodes_in_col = nodes_in_col
        self.nodes_in_row = nodes_in_row
        self.empty_rows = empty_rows
        self.boundary_conditions = boundary_conditions

        self.edges_in_row = 3*self.nodes_in_row - 1
        self.edges = np.array()

    def neighbor(self, i: int, j: int, k: int) -> Tuple[int, int]:
        neighbors = util.Neighbors()
        if j % 3 == 0:          # We are a y-edge
            if j == 0:
                j += neighbors.LEFT[k]
            if j == self.edges_in_row-3:
                j += neighbors.RIGHT[k]
            i += neighbors.Y[0][k]
            j += neighbors.Y[1][k]
            return i, j
        if j % 3 == 1:          # We are a z-edge
            if i == 0:
                i += neighbors.BOTTOM[k]
            if i == self.nodes_in_row-1:
                i += neighbors.TOP[k]
            if j == 1:
                j += neighbors.LEFT[k]
            if j == self.edges_in_row-2:
                j += neighbors.RIGHT[k]
            i += neighbors.Z[0][k]
            j += neighbors.Z[1][k]
            return i, j
        if j % 3 == 2:          # We are an x-edge
            if i == 0:
                i += neighbors.BOTTOM[k]
            if i == self.nodes_in_row-1:
                i += neighbors.TOP[k]
            i += neighbors.X[0][k]
            j += neighbors.X[1][k]
            return i, j
        pass
