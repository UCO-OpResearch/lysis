from enum import Enum, unique

import numpy as np





class FiberGrid(object):
    def __init__(self, height, width, empty_rows, boundary_conditions):
        self.height = height
        self.width = width
        self.empty_rows = empty_rows
        self.boundary_conditions = boundary_conditions



