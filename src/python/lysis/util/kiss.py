"""Minimal skeleton for using the Marsaglia KISS generator

Wrapper for C module.

Example Usage:

* Initialize
    >>> kiss = KissGenerator()
    >>> kiss.set_seed(123)

* Generate random integer 0..(2^32)-1
    >>> kiss.kiss32()

* Generate random U[0,1]
    >>> kiss.urcw1()
"""

import ctypes
import os
from typing import Tuple


class KissGenerator(object):
    """A Pseudo-random Number Generator.

    Wrapper class for kiss.so C module
    """
    def __init__(self, seed: int = None, state: Tuple[int, int, int, int] = None):
        """Creates a new KISS random number generator object

        If no seed is given, one is generated from the system clock.

        Args:
            seed: The seed for the RNG
            state: The state of the RNG
        """
        # Determine the current path and find the kiss.so library
        path = os.path.dirname(__file__)
        lib_path = os.path.join(path, '..', '..', '..', '..', 'lib')
        kiss_file = 'kiss.so'

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
        # It accepts a (pointer to a) 4-element array of unsigned integers into which it writes the state
        self.__get_kiss32.argtypes = [self.state_type]

        if state is not None:
            # If a state was given, then set it
            self.set_state(state)
        elif seed is not None:
            # If there was no state given, but a seed instead, use the seed.
            self.set_seed(seed)
        else:
            # If there was no state or seed given, generate a seed from the system clock and use that.
            self.set_seed(self.mscw())
        
    def set_state(self, state: Tuple[int, int, int, int]):
        """Sets the state of the Random Number Generator.

        Args:
            state: A tuple of four integers.

                Note: While Python will accept any 64-bit integer (-9,223,372,036,854,775,806 through
                9,223,372,036,854,775,807), these will be converted to unsigned 32-bit integers (0 through
                4,294,967,295) when passed to the underlying C code with unpredictable results.
        """
        # Unpack the state tuple into a c-type array
        c_state = self.state_type(state[0], state[1], state[2], state[3])
        # Pass the state to the C generator
        self.__set_kiss32(c_state)

    def set_seed(self, seed: int):
        """Sets the seed for the Random Number Generator.

        Args:
            seed: The seed.

                Note: While Python will accept any 64-bit integer (-9,223,372,036,854,775,806 through
                9,223,372,036,854,775,807), this will be converted to an unsigned 32-bit integer (0 through
                4,294,967,295) when passed to the underlying C code with unpredictable results.
        """
        # Get the current state (since the seed is the fourth element of the state)
        state = self.get_state()
        # Replace the fourth element of the state with the new seed
        state = state[:3] + (seed,)
        # Set the new state
        self.set_state(state)

    def get_state(self) -> Tuple[int, int, int, int]:
        """Returns the current state of the Random Number Generator."""
        # Define a C array to hold the state
        c_state = self.state_type(0, 0, 0, 0)
        # Send the C array to the module
        # Note that this is a C-style "pass by reference". That is, the function does not "return" anything, but the
        # memory that the "c_state" pointer points to will be modified by the function call.
        self.__get_kiss32(c_state)
        # Unpack the C array and convert to a tuple
        return tuple(int(i) for i in c_state)

    def urcw1(self) -> float:
        """Returns the next double-precision random number from the pseudo-random number stream. This is distributed
        uniformly from zero through one."""
        # This method will be overwritten by the one from the C library when the class is initialized.
        pass

    def mscw(self) -> int:
        """Returns a pseudo-random 32-bit unsigned integer (i.e., [0..(2^32)-1]).
        This is based on the system clock and NOT the current seed or state."""
        # This method will be overwritten by the one from the C library when the class is initialized.
        pass

    def kiss32(self) -> int:
        """Returns the next 32-bit unsigned integer from the pseudo-random number stream. This is distributed uniformly
        from 0 through (2^32)-1."""
        # This method will be overwritten by the one from the C library when the class is initialized.
        pass
