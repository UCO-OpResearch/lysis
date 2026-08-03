import numpy as np
from ...config.run import Run

class BindingTimeFactory:
    """Factory for efficient batch generation of molecule binding times.

    This class pre-generates binding times in large batches to amortize the
    cost of random number generation. It maintains an internal buffer that is
    automatically refilled when exhausted, providing better performance than
    generating binding times one at a time or in small batches.

    Binding times are calculated as:
        binding_time = -ln(X) / (binding_rate * binding_sites) - time_step/2
    where X ~ Uniform(0,1).

    Note: The current time is NOT included in generated binding times and must
    be added by the caller.

    :param run: The Run object containing simulation parameters
    :type run: Run
    :param rng: Random number generator shared with the main simulation to
        prevent statistical overlap
    :type rng: np.random.Generator or KissRandomGenerator
    """

    def __init__(self, run: Run, rng: np.random.Generator):
        # Buffer size: 1M minimum, or 10x molecule count (to amortize refill costs)
        self.period = max(1_000_000, 10 * run.macro_params.total_molecules)
        # The run this object is a part of
        self.run = run
        # The random number generator to be used
        self.rng = rng
        # Initialize the list. It starts empty and will be filled the first
        # time `next()` is called.
        self.list = np.empty((self.period,), dtype=np.double)
        # This points to the next binding time in the list to be served
        # We initialize it here at the end of the list as the list is
        # currently empty
        self.pointer = self.period

    def fill_list(
        self,
        random_numbers: np.ndarray = None,
    ) -> np.ndarray | None:
        """Generate binding times from random numbers.

        This method operates in two modes:

        1. Internal mode (random_numbers=None): Refills the internal buffer with
            new binding times. This is the standard mode that provides amortized
            performance benefits. Returns None.

        2. External mode (random_numbers provided): Calculates binding times for
            the provided random numbers without affecting the internal buffer.
            Returns the calculated binding times. Does NOT benefit from amortization.

        Binding times are calculated as:
            binding_time = -ln(X) / (binding_rate * binding_sites) - time_step/2
        where X values come from either the RNG or the provided array.

        Note: Generated binding times do NOT include the current time. The caller
        must add current_time to get absolute binding times.

        :param random_numbers: Optional array of random numbers in [0,1) to convert
            to binding times. If None, generates new random numbers and refills
            the internal buffer.
        :type random_numbers: np.ndarray, optional
        :return: Array of binding times if random_numbers was provided, otherwise None
        :rtype: np.ndarray | None
        """
        # Check to see if we were passed any random numbers
        if random_numbers is None:
            # If not, we need to generate random numbers
            binding_time_list = self.rng.random(len(self.list))
            # Set a flag so that we know to use the internal list
            return_list = False
        else:
            # If we were passed random numbers, we use those
            binding_time_list = random_numbers
            # And set a flag so that we know to return the times we generate
            return_list = True

        # Binding time is calculated as:
        # (current_time) - ln(X)/(binding_rate * binding_sites) - time_step/2
        # where X ~ U(0,1)
        # In general, we don't know the current time when the binding time
        # is generated, so we depend on the calling code to add that.

        # Calculate binding time = -ln(X) / (rate * sites) - timestep/2
        # All operations done in-place for memory efficiency
        binding_time_list = np.log(binding_time_list, out=binding_time_list)
        denominator = (
            self.run.macro_params.micro_params.bind_rate_tPA
            * self.run.macro_params.micro_params.binding_sites
        ).magnitude
        binding_time_list = np.divide(
            binding_time_list, denominator, out=binding_time_list
        )
        # Use subtract instead of negating then adding (equivalent but cleaner)
        binding_time_list = np.subtract(
            -self.run.macro_params.time_step.magnitude / 2,
            binding_time_list,
            out=binding_time_list,
        )
        # Either return the binding times we calculated or store them in
        # the internal list.
        if return_list:
            return binding_time_list
        else:
            self.list = binding_time_list

    def next(self, count: int = 1):
        """Provide the next batch of binding times from the internal buffer.

        Serves binding times from the internal buffer, automatically refilling it
        when necessary. The running time is inconsistent (fast most of the time,
        slow when refilling), but amortized performance is excellent.

        This is the preferred method for obtaining binding times during simulation
        as it provides optimal performance through batching.

        :param count: Number of binding times requested
        :type count: int
        :return: Array of binding times (excluding current time, which must be added by caller)
        :rtype: np.ndarray
        :raises: May raise memory errors if count exceeds available memory
        """
        # Allocate the array for the output
        out = np.empty(count, dtype=np.double)
        # If there are enough binding times left in the internal array to
        # provide what is being requested
        if count + self.pointer <= self.period:
            # Fill the output array from the internal array
            out = self.list[self.pointer : self.pointer + count]
            # Update the pointer
            self.pointer += count
        # Otherwise, if there are not enough stored binding times
        else:
            # Not enough binding times left in buffer - need to span across refill
            remaining = self.period - self.pointer
            # Copy remaining times from old buffer
            out[:remaining] = self.list[self.pointer :]
            # Refill the internal buffer
            self.fill_list()
            # Copy additional needed times from newly filled buffer
            needed_from_new = count - remaining
            out[remaining:] = self.list[:needed_from_new]
            # Update pointer to position after times we just served
            self.pointer = needed_from_new
        return out