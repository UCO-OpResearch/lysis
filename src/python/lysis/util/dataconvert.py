from dataclasses import asdict
from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Mapping, Union, Callable


import numpy as np
import h5py

from .constants import CONST
from .dataspec import (
    DataCollectionSpec,
    DataCollectionType,
    DataSetSpec,
    DataSetType,
    DataSpec,
    dataspec,
    parse_shape,
    check_dataset_spec,
)
from .edge_grid import generate_fortran_neighborhood_structure

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


data_converters: dict[
    tuple[DataSpec, DataSpec], dict[str, Callable[[DataCollectionType], DataSetType]]
] = {
    (dataspec["v2.0.0"], dataspec["v1.99.0"]): {
        "micro_log": lambda data: data["micro_log"],
        "firstPLi": lambda data: data["pli_first_time"],
    }
}


def generate_macroscale_in(in_data: DataSetType) -> DataSetType:
    """
    This method generates the input data for the Fortran macroscale model.
    *NOTE: Even though this method should only be used with the Fortran macroscale model,
    the code itself uses the HDF5 data specification.
    The output MUST be converted to Fortran specification before writing to disk.*

    :param in_data: _description_
    :type in_data: DataSetType
    :return: _description_
    :rtype: DataSetType
    """
    out_data = in_data.copy()
    # Get the number of microscale runs and set the dimensions of the bins so that we get 100 bins
    bin_size = in_data["pli_first_time"].size // 100
    # bin_edge_proportions is the CDF of the tPA leaving time distribution.
    # This is really just a list of edgepoints from the bins for tPA leaving time
    # These bins are evenly distributed along the interval [0, 1]
    out_data["bin_edge_proportions"] = np.append(np.arange(0, 1, 0.01), [1.0])

    # The remaining data will be arranged into 100 bins
    # according to the time tPA left the simulation.
    # Get the sorted ordering of the tPA leaving times
    indices = in_data["tpa_leaving_time"].argsort()
    # Find the tPA leaving times for the edges of each bin.
    out_data["bin_edge_tpa_leaving_time"] = np.append(
        [0], in_data["tpa_leaving_time"][indices[bin_size - 1 :: bin_size]]
    )

    # Identify which simulations had the fiber fully degraded
    lysis_complete = in_data["fiber_degraded"]
    # Read in the fiber degradation times for all simulations
    lysis_time = in_data["sim_final_time"]
    # If full degradation did NOT occur,
    # this matrix currently contains the ending time of the simulation.
    # Replace these times with an 'infinity' marker of 6,000 seconds
    lysis_time[~lysis_complete] = 6_000
    # Rearrange the matrix so that each row contains the lysis times for simulations
    # corresponding to the matching bin in the ``bin_edge_proportions`` vector.
    # Then sort the rows (bins) individually by lysis time.
    # Finally, transpose the matrix so that the bins are arranged in columns.
    out_data["binned_fiber_degrade_time"] = np.stack(
        [
            np.sort(lysis_time[indices[i * bin_size : (i + 1) * bin_size]])
            for i in range(100)
        ]
    ).T

    # Find the location of the first '6000' entry in each column of the ``binned_fiber_degrade_time`` matrix
    # Then convert to 1-indexing.
    out_data["binned_fiber_degraded"] = out_data["binned_fiber_degrade_time"].argmax(
        axis=0
    )

    # Generate the list of neighborhoods for the fortran code.
    out_data["edge_grid_neighbors"] = generate_fortran_neighborhood_structure(
        in_data["params"]["macro_params"]["rows"],
        in_data["params"]["macro_params"]["cols"],
    )


def convert_data(
    input_data: DataCollectionType,
    input_set_spec: DataSetSpec,
    output_set_spec: DataSetSpec,
) -> DataCollectionType:
    out_data = {}
    for collection in output_set_spec:
        for dataset_needed in output_set_spec[collection].data.keys():
            try:
                out = data_converters[input_set_spec, output_set_spec][dataset_needed](
                    input_data
                )
            except KeyError:
                continue
            else:
                out_data[dataset_needed] = out
    return out_data
