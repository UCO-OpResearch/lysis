import functools
import warnings

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
    dataspec,
    parse_shape,
    check_dataset_spec,
    tags,
)
from .edge_grid import (
    generate_fortran_neighborhood_structure,
    from_fortran_edge_index_array,
)

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


def safe_np_int_conversion(int_array, dtype=np.uint8, copy=True):
    """
    A few functions require arrays of a certain type (e.g. np.int32, np.uint8).
    To allow functions to accept standard numpy integer arrays (usually of
    dtype=np.int64) we cast but check bounds to avoid wrap-around
    conversion errors (numpy doesn't seem to provide this functionality)

    Shamelessly stolen from https://stackoverflow.com/questions/56684893/numpy-cast-from-signed-to-unsigned-int-with-same-kind
    Modified to allow boolean
    """
    int_array = np.array(int_array)
    if int_array.size == 0:
        return int_array.astype(dtype, copy=copy)  # Allow empty arrays of any type
    try:
        return int_array.astype(dtype, casting="safe", copy=copy)
    except TypeError:
        bounds = np.iinfo(dtype)
        if np.all(int_array >= bounds.min) and np.all(int_array <= bounds.max):
            if int_array.dtype.kind == "i" and np.dtype(dtype).kind == "u":
                # Allow casting from int to unsigned int, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from, e.g., a float.
                casting = "same_kind"
            return int_array.astype(dtype, casting=casting, copy=copy)
        else:
            raise OverflowError("Cannot convert safely to {} type".format(dtype))


def safe_np_bool_conversion(int_array, copy=True):
    """
    A few functions require arrays of np.bool.
    To allow functions to accept standard numpy integer arrays (usually of
    dtype=np.int64) we cast but check bounds to avoid wrap-around
    conversion errors (numpy doesn't seem to provide this functionality)

    Shamelessly stolen from https://stackoverflow.com/questions/56684893/numpy-cast-from-signed-to-unsigned-int-with-same-kind
    Modified to allow boolean
    """
    int_array = np.array(int_array)
    if int_array.size == 0:
        return int_array.astype(np.bool, copy=copy)  # Allow empty arrays of any type
    try:
        return int_array.astype(np.bool, casting="safe", copy=copy)
    except TypeError:
        if np.all(int_array >= 0) and np.all(int_array <= 1):
            if int_array.dtype.kind == "i":
                # Allow casting from int to bool, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from, e.g., a float.
                casting = "same_kind"
            return int_array.astype(np.bool, casting=casting, copy=copy)
        else:
            raise OverflowError("Cannot convert safely to np.bool type")


def generate_macroscale_in(in_data: DataCollectionType) -> DataCollectionType:
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

    # TODO: Add code to check that `in_data` meets the "current" data specification
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

    # Calculate
    out_data["params"]["macro_params"]["forced_unbind"] = np.count_nonzero(
        in_data["tPA_forced_unbind"]
    ) / (
        np.count_nonzero(in_data["tPA_forced_unbind"])
        + np.count_nonzero(in_data["tPA_kinetic_unbind"])
    )

    return out_data


def convert_fiber_degrade_time(
    input_data: DataCollectionType, input_spec: str, output_spec: str
) -> DataSetType:
    """
    TODO: Function from cell 9 of H5-File-Builder.ipynb
    _summary_

    :param data: a dict of data in fortran format
    :type data: DataCollectionType
    :return: A numpy array, or a list of numpy arrays in `output_spec` format
    """
    # TODO: Check if simulations_combined = True
    output_data = []
    for data in input_data["f_deg_list"]:
        output_data.append(
            np.empty(
                data.shape,
                dtype=dataspec["v2.0.0"]["macroscale_out"]
                .data["fiber_degrade_time"]
                .dtype,
            )
        )
        output_data[-1][["Simulation Time Elapsed", "Fiber New Degrade Time"]] = data[
            ["Simulation Time Elapsed", "Fiber New Degrade Time"]
        ]
        locations = from_fortran_edge_index_array(
            data["Grid Location Index"] - 1,
            input_data["params"]["macro_params"]["rows"],
            input_data["params"]["macro_params"]["cols"],
        )
        output_data[-1]["Grid Location Row"] = locations[:, 0]
        output_data[-1]["Grid Location Rank"] = locations[:, 1]
    return output_data


# TODO Do the same thing with the tpa_bind_events from cell 10 of H5-File-Builder.ipynb


def _not_implemented(dataset_name: str, input_spec: str, output_spec: str):
    """Create a converter stub that raises NotImplementedError with a clear message.

    :param dataset_name: Name of the dataset being converted
    :type dataset_name: str
    :param input_spec: Input specification version
    :type input_spec: str
    :param output_spec: Output specification version
    :type output_spec: str
    :return: Function that raises NotImplementedError when called
    :rtype: Callable
    """
    def _raise_error(data):
        raise NotImplementedError(
            f"Conversion for dataset '{dataset_name}' from {input_spec} to {output_spec} "
            f"is not yet implemented."
        )
    return _raise_error


data_converters: dict[
    tuple[str, str], dict[str, Callable[[DataCollectionType], DataSetType]]
] = {
    ("v2.0.0", "v1.99.0"): {
        # Microscale output (implemented)
        "micro_log": lambda data: data["micro_log"],
        "firstPLi": lambda data: data["pli_first_time"],
        "lasttPA": lambda data: data["tpa_final_num"],
        "lyscomplete": lambda data: data["fiber_degraded"],
        "lysis": lambda data: data["sim_final_time"],
        "PLi": lambda data: data["pli_generated_num"],
        "tPA_time": lambda data: data["tpa_leaving_time"],
        "tPAPLiunbd": lambda data: data["tpa_unbound_by_pli"],
        "tPAunbind": lambda data: data["tpa_unbound_kinetic"],
        # Macroscale input (not implemented)
        "tPAleave": _not_implemented("tPAleave", "v2.0.0", "v1.99.0"),
        "tsectPA": _not_implemented("tsectPA", "v2.0.0", "v1.99.0"),
        "lysismat": _not_implemented("lysismat", "v2.0.0", "v1.99.0"),
        "lenlysisvect": _not_implemented("lenlysisvect", "v2.0.0", "v1.99.0"),
        "neighbors": _not_implemented("neighbors", "v2.0.0", "v1.99.0"),
        # Macroscale output (not implemented)
        "macro_log": lambda data: data["macro_log"],
        "Nsave": _not_implemented("Nsave", "v2.0.0", "v1.99.0"),
        "tsave": _not_implemented("tsave", "v2.0.0", "v1.99.0"),
        "f_deg_list": _not_implemented("f_deg_list", "v2.0.0", "v1.99.0"),
        "m_bind_t": _not_implemented("m_bind_t", "v2.0.0", "v1.99.0"),
        "m_loc": _not_implemented("m_loc", "v2.0.0", "v1.99.0"),
        "m_bound": _not_implemented("m_bound", "v2.0.0", "v1.99.0"),
        "mfpt": _not_implemented("mfpt", "v2.0.0", "v1.99.0"),
    },
    ("v1.99.0", "v2.0.0"): {
        # Microscale output (implemented)
        "micro_log": lambda data: data["micro_log"],
        "pli_first_time": lambda data: data["firstPLi"],
        "tpa_final_num": lambda data: data["lasttPA"],
        "fiber_degraded": lambda data: data["lyscomplete"],
        "sim_final_time": lambda data: data["lysis"],
        "pli_generated_num": lambda data: data["PLi"],
        "tpa_leaving_time": lambda data: data["tPA_time"],
        "tpa_unbound_by_pli": lambda data: data["tPAPLiunbd"],
        "tpa_unbound_kinetic": lambda data: data["tPAunbind"],
        # Macroscale input (not implemented)
        "bin_edge_proportions": _not_implemented("bin_edge_proportions", "v1.99.0", "v2.0.0"),
        "bin_edge_tpa_leaving_time": _not_implemented("bin_edge_tpa_leaving_time", "v1.99.0", "v2.0.0"),
        "binned_fiber_degrade_time": _not_implemented("binned_fiber_degrade_time", "v1.99.0", "v2.0.0"),
        "binned_fiber_degraded": _not_implemented("binned_fiber_degraded", "v1.99.0", "v2.0.0"),
        "edge_grid_neighbors": _not_implemented("edge_grid_neighbors", "v1.99.0", "v2.0.0"),
        # Macroscale output (partially implemented)
        "macro_log": lambda data: data["macro_log"],
        "snapshot_time": _not_implemented("snapshot_time", "v1.99.0", "v2.0.0"),
        "fiber_degrade_time": functools.partial(
            convert_fiber_degrade_time, input_spec="v1.99.0", output_spec="v2.0.0"
        ),
        "tpa_bind_events": _not_implemented("tpa_bind_events", "v1.99.0", "v2.0.0"),
        "tpa_location_snapshot": _not_implemented("tpa_location_snapshot", "v1.99.0", "v2.0.0"),
        "tpa_transit_time": _not_implemented("tpa_transit_time", "v1.99.0", "v2.0.0"),
    },
}


def convert_data(
    input_data: DataCollectionType,
    input_set_spec: str,
    output_set_spec: str,
) -> DataCollectionType:
    while input_set_spec in tags:
        input_set_spec = tags[input_set_spec]
    while output_set_spec in tags:
        output_set_spec = tags[output_set_spec]
    # TODO: Check that `input_set_spec` and `output_set_spec` exist
    # TODO: Add code to check that `input_data` meets the specifications of `input_set_spec`.
    out_data = {}
    out_data["params"] = input_data["params"]
    for collection in dataspec[output_set_spec].values():
        for dataset_needed in collection.data.keys():
            try:
                out = data_converters[input_set_spec, output_set_spec][dataset_needed](
                    input_data
                )
            except KeyError as e:
                warnings.warn(f"Missing data for {dataset_needed}")
                continue
            if collection.simulations_combined:
                out = [out]
            out_data[dataset_needed] = []
            for data_table in out:
                try:
                    output_table = data_table.astype(
                        collection.data[dataset_needed].dtype,
                        casting="same_kind",
                    )
                except TypeError as e:
                    dt = np.dtype(collection.data[dataset_needed].dtype)
                    match dt.kind:
                        case "u":
                            output_table = safe_np_int_conversion(data_table, dtype=dt)
                        case "b":
                            output_table = safe_np_bool_conversion(data_table)
                        case "f":
                            # TODO: Create a function that does the same thing as the safe_np_int_conversion, but for floats.
                            raise NotImplementedError("Not implemented yet")
                        case _:
                            raise e
                out_data[dataset_needed].append(output_table)
            if collection.simulations_combined:
                out_data[dataset_needed] = out_data[dataset_needed][0]
    return out_data
