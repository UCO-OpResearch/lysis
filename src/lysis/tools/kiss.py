"""Python wrapper for the Marsaglia KISS random number generator.

This module provides a Python interface to the KISS (Keep It Simple Stupid)
random number generator implemented in C. KISS is a combination of multiple
simple generators that produces high-quality pseudo-random numbers with a
very long period.

The KISS generator combines:
- A 3-shift-register generator (xorshift)
- A multiply-with-carry generator
- A congruential generator

This implementation is used primarily for validation against legacy Fortran
code that uses the same KISS algorithm.

Example Usage:

Initialize with a seed::

    >>> kiss = KissRandomGenerator()
    >>> kiss.seed(123)

Generate random integer in range [0, 2^32-1]::

    >>> kiss.kiss32()

Generate random float in range [0, 1]::

    >>> kiss.random()

Generate array of random floats::

    >>> kiss.random(size=100)
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
    """Pseudo-random number generator using the KISS algorithm.

    This class wraps the KISS (Keep It Simple Stupid) random number generator
    implemented in C (kiss.so shared library). It provides a numpy-compatible
    interface for generating random numbers using the same algorithm as the
    original Fortran implementation.

    The KISS generator has excellent statistical properties and a period of
    approximately 2^123, making it suitable for Monte Carlo simulations.

    The interface mimics numpy.random.Generator to allow drop-in replacement
    in code that needs to exactly replicate the Fortran KISS sequence.

    Attributes:
        state_type: ctypes array type for the 4-element generator state
        urcw1: C function for generating single U(0,1) random numbers
        vurcw1: C function for generating arrays of U(0,1) random numbers
        mscw: C function for generating time-based seeds
        kiss32: C function for generating 32-bit unsigned integers

    Note:
        The underlying C library (kiss.so) must be compiled and available in
        the lib/ directory relative to this module.
    """

    def __init__(self, seed: int = None):
        """Initialize a new KISS random number generator.

        Loads the kiss.so C library, imports all necessary C functions, sets up
        the ctypes interfaces, and initializes the generator state with the
        provided seed (or a time-based seed if none is provided).

        :param seed: Initial seed for the RNG. If None, generates a seed from
            the system clock. Must be in range [0, 2^32-1] (values outside this
            range will be truncated when passed to C).
        :type seed: int, optional
        :raises OSError: If the kiss.so library cannot be found or loaded
        :raises ctypes.CException: If C function signatures don't match expected types
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

        self.vurcw1 = my_kiss.c_vurcw2_
        self.vurcw1.argtypes = [
            np.ctypeslib.ndpointer(np.float_, flags="C_CONTIGUOUS"),
            ctypes.c_int,
        ]

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
        """Set the internal state of the random number generator.

        This allows restoration of a previously saved generator state, enabling
        reproducible random sequences. The state consists of four 32-bit unsigned
        integers that encode the current position in the KISS sequence.

        :param state: A tuple of four integers representing the generator state.
            Each integer should be in the range [0, 2^32-1]. Values outside this
            range will be truncated to 32 bits with unpredictable results.
        :type state: Tuple[int, int, int, int]

        Warning:
            Python accepts 64-bit integers (-9,223,372,036,854,775,808 through
            9,223,372,036,854,775,807), but these will be converted to unsigned
            32-bit integers (0 through 4,294,967,295) when passed to the C code.
            Values outside the 32-bit range will produce undefined behavior.

        Example:
            >>> kiss = KissRandomGenerator(123)
            >>> state = kiss.getstate()
            >>> kiss.random()  # Generate some numbers
            >>> kiss.setstate(state)  # Restore previous state
        """
        # Unpack the state tuple into a c-type array
        c_state = self.state_type(state[0], state[1], state[2], state[3])
        # Pass the state to the C generator
        self.__set_kiss32(c_state)

    def seed(self, seed: int = None):
        """Set the seed for the random number generator.

        The seed determines the starting point of the random sequence. The same
        seed will always produce the same sequence of random numbers, enabling
        reproducible simulations.

        Internally, the seed becomes the fourth element of the generator's state
        vector. The other three state elements are preserved from the current state.

        :param seed: The seed value. If None, generates a time-based seed using
            the system clock. Should be in range [0, 2^32-1].
        :type seed: int, optional

        Warning:
            Python accepts 64-bit integers (-9,223,372,036,854,775,808 through
            9,223,372,036,854,775,807), but these will be converted to unsigned
            32-bit integers (0 through 4,294,967,295) when passed to the C code.
            Values outside the 32-bit range will produce undefined behavior.

        Example:
            >>> kiss1 = KissRandomGenerator(42)
            >>> kiss2 = KissRandomGenerator(42)
            >>> kiss1.random() == kiss2.random()  # Same seed, same sequence
            True
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
        """Get the current internal state of the random number generator.

        Returns the four 32-bit unsigned integers that represent the current
        position in the KISS sequence. This state can be saved and later restored
        using setstate() to resume generation from the same point.

        :return: A tuple of four integers representing the generator state.
            Each integer is in the range [0, 2^32-1].
        :rtype: Tuple[int, int, int, int]

        Example:
            >>> kiss = KissRandomGenerator(123)
            >>> state = kiss.getstate()
            >>> print(state)
            (123, 234567890, 345678901, 456789012)  # Example values
        """
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
        """Generate random float(s) uniformly distributed in [0, 1).

        If size is None, returns a single float. If size is specified, returns
        a numpy array of the given size filled with random floats.

        This method mimics numpy.random.Generator.random() for compatibility.

        :param size: If None, return a single float. If an integer, return a
            1D numpy array of that length.
        :type size: int, optional
        :return: A single random float in [0, 1) if size is None, otherwise a
            numpy array of random floats.
        :rtype: float or np.ndarray

        Example:
            >>> kiss = KissRandomGenerator(42)
            >>> kiss.random()  # Single value
            0.123456789
            >>> kiss.random(5)  # Array of 5 values
            array([0.234, 0.567, 0.890, 0.123, 0.456])
        """
        if size is None:
            return self.urcw1()
        else:
            out = np.empty(size, dtype=np.float_)
            self.vurcw1(out, size)
            # for i in range(size):
            #     out[i] = self.urcw1()
            return out

    def integers(
        self, bottom: int, top: int = None, size: int = None
    ) -> int | np.ndarray:
        """Generate random integer(s) in a specified range.

        Returns random integers from the half-open interval [bottom, top). If
        top is not specified, returns integers from [0, bottom).

        This method mimics numpy.random.Generator.integers() for compatibility.

        :param bottom: If top is None, this is the exclusive upper bound and the
            lower bound is 0. If top is specified, this is the inclusive lower bound.
        :type bottom: int
        :param top: Exclusive upper bound. If None, bottom becomes the upper bound
            and 0 becomes the lower bound.
        :type top: int, optional
        :param size: If None, return a single integer. If an integer, return a
            1D numpy array of that length.
        :type size: int, optional
        :return: A single random integer if size is None, otherwise a numpy array
            of random integers in the range [bottom, top).
        :rtype: int or np.ndarray

        Example:
            >>> kiss = KissRandomGenerator(42)
            >>> kiss.integers(10)  # Random int in [0, 10)
            7
            >>> kiss.integers(5, 15)  # Random int in [5, 15)
            12
            >>> kiss.integers(0, 10, size=5)  # Array of 5 random ints
            array([3, 7, 2, 9, 1])
        """
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
        """Generate a time-based pseudo-random 32-bit unsigned integer.

        This function uses the system clock to generate a random seed value,
        making it useful for initialization when reproducibility is not required.

        Important: This function is INDEPENDENT of the generator's current state
        and seed. It always uses the system clock, not the KISS sequence.

        :return: A pseudo-random integer in the range [0, 2^32-1] based on the
            current system time.
        :rtype: int

        Note:
            This method stub is overwritten by the C library function during
            __init__. The implementation here is never actually called.
        """
        # This method will be overwritten by the one from the C library when
        # the class is initialized.
        pass

    def kiss32(self) -> int:
        """Generate the next 32-bit unsigned integer from the KISS sequence.

        This is the core KISS generator function that produces uniformly
        distributed integers in the full 32-bit range. The random() method uses
        this internally and converts the output to floats in [0, 1).

        :return: The next random integer in the range [0, 2^32-1], uniformly
            distributed.
        :rtype: int

        Note:
            This method stub is overwritten by the C library function during
            __init__. The implementation here is never actually called.

        Example:
            >>> kiss = KissRandomGenerator(42)
            >>> kiss.kiss32()
            3141592653  # Example value
        """
        # This method will be overwritten by the one from the C library when
        # the class is initialized.
        pass
