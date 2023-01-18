"""Minimal skeleton for using the Marsaglia KISS generator

Wrapper for C module.

Example Usage:

* Initialize
    >>> kiss = KissRandomGenerator()
    >>> kiss.seed(123)

* Generate random integer 0..(2^32)-1
    >>> kiss.kiss32()

* Generate random U[0,1]
    >>> kiss.random()
"""

import ctypes
import os
from typing import Tuple

import numpy as np

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class KissRandomGenerator:
    """A Pseudo-random Number Generator.

    Wrapper class for kiss.so C module
    """

    def __init__(self, seed: int = None):
        """Creates a new KISS random number generator object

        If no seed is given, one is generated from the system clock.

        Args:
            seed: The seed for the RNG
        """
        # Determine the current path and find the kiss.so library
        path = os.path.dirname(__file__)
        lib_path = os.path.join(path, "..", "..", "..", "..", "lib")
        kiss_file = "kiss.so"

        # Import the C library
        my_kiss = ctypes.CDLL(os.path.join(lib_path, kiss_file))

        # Defile the datatype for the RNG state
        self.state_type = ctypes.c_uint * 4

        # Import the U(0,1) generator function
        self.urcw1 = my_kiss.urcw1_
        # It returns a double-precision (64-bit) float (equivalent to Python float)
        self.urcw1.restype = ctypes.c_double

        # Import the random seed function
        self.mscw = my_kiss.mscw_
        # It returns an unsigned 32-bit integer
        self.mscw.restype = ctypes.c_uint32

        # Import the random integer generator
        self.kiss32 = my_kiss.kiss32_
        # It returns an unsigned 32-bit integer
        self.kiss32.restype = ctypes.c_uint32

        # Import the method to set the state
        self.__set_kiss32 = my_kiss.set_kiss32_
        # It accepts a 4-element array of unsigned integers
        self.__set_kiss32.argtypes = [self.state_type]

        # Import the method to get the state
        self.__get_kiss32 = my_kiss.get_kiss32_
        # It accepts a (pointer to a) 4-element array of unsigned integers into
        # which it writes the state
        self.__get_kiss32.argtypes = [self.state_type]

        # Set the seed with what was given
        self.seed(seed)

    def setstate(self, state: Tuple[int, int, int, int]):
        """Sets the state of the Random Number Generator.

        Args:
            state: A tuple of four integers.

                Note: While Python will accept any 64-bit integer
                (-9,223,372,036,854,775,806 through 9,223,372,036,854,775,807),
                these will be converted to unsigned 32-bit integers (0 through
                4,294,967,295) when passed to the underlying C code with
                unpredictable results.
        """
        # Unpack the state tuple into a c-type array
        c_state = self.state_type(state[0], state[1], state[2], state[3])
        # Pass the state to the C generator
        self.__set_kiss32(c_state)

    def seed(self, seed: int = None):
        """Sets the seed for the Random Number Generator.

        Args:
            seed: The seed.

                Note: While Python will accept any 64-bit integer
                (-9,223,372,036,854,775,806 through 9,223,372,036,854,775,807),
                this will be converted to an unsigned 32-bit integer (0 through
                4,294,967,295) when passed to the underlying C code with
                unpredictable results.
        """
        # If no seed was given, generate one from the system clock
        if seed is None:
            seed = self.mscw()
        # Get the current state (since the seed is the fourth element of the state)
        state = self.getstate()
        # Replace the fourth element of the state with the new seed
        state = (int(state[0]), int(state[1]), int(state[2]), seed)
        # Set the new state
        self.setstate(state)

    def getstate(self) -> Tuple[int, int, int, int]:
        """Returns the current state of the Random Number Generator."""
        # Define a C array to hold the state
        c_state = self.state_type(0, 0, 0, 0)
        # Send the C array to the module
        # Note that this is a C-style "pass by reference". That is, the function
        # does not "return" anything, but the memory that the "c_state" pointer
        # points to will be modified by the function call.
        self.__get_kiss32(c_state)
        # Unpack the C array and convert to a tuple
        return c_state[0], c_state[1], c_state[2], c_state[3]

    def random(self, size: int = None) -> float | np.ndarray:
        """Returns the next double-precision random number from the
        pseudo-random number stream. This is distributed uniformly from zero
        through one."""
        if size is None:
            return self.urcw1()
        else:
            out = np.empty((size,), dtype=int)
            for i in range(size):
                out[i] = self.urcw1()
            return out

    def integers(
        self, bottom: int, top: int = None, size: int = None
    ) -> int | np.ndarray:
        if top is None:
            top = bottom
            bottom = 0
        if size is None:
            return bottom + int((top - bottom) * self.urcw1())
        else:
            out = np.empty((size,), dtype=int)
            for i in range(size):
                out[i] = bottom + int((top - bottom) * self.urcw1())
            return out

    def mscw(self) -> int:
        """Returns a pseudo-random 32-bit unsigned integer
        (i.e., [0..(2^32)-1]).
        This is based on the system clock and NOT the current seed or state."""
        # This method will be overwritten by the one from the C library when
        # the class is initialized.
        pass

    def kiss32(self) -> int:
        """Returns the next 32-bit unsigned integer from the pseudo-random
        number stream. This is distributed uniformly from 0 through (2^32)-1."""
        # This method will be overwritten by the one from the C library when
        # the class is initialized.
        pass
