import os
import secrets
import time

from math import log, floor
from typing import Any, AnyStr, Mapping


tokens = "23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"
_last_v8_timestamp = None


def uuid8code():
    num = uuid8()
    length = floor(128 * log(2) / log(len(tokens))) + 1
    result = ""
    while num > 0:
        num, remainder = divmod(num, len(tokens))
        result += tokens[remainder]
    for i in range(length - len(result)):
        result += "0"
    result = result[::-1]
    result = "-".join(
        [result[:4], result[4:8], result[8:12], result[12:16], result[16:]]
    )
    return result


def uuid8():
    r"""UUID version 8 features a time-ordered value field derived from the
    widely implemented and well known Unix Epoch timestamp source, the
    number of nanoseconds seconds since midnight 1 Jan 1970 UTC, leap
    seconds excluded.

    Copied from https://github.com/oittaa/uuid6-python under MIT License
    """

    global _last_v8_timestamp

    nanoseconds = time.time_ns()
    if _last_v8_timestamp is not None and nanoseconds <= _last_v8_timestamp:
        nanoseconds = _last_v8_timestamp + 1
    _last_v8_timestamp = nanoseconds
    timestamp_ms, timestamp_ns = divmod(nanoseconds, 10**6)
    subsec = timestamp_ns * 2**20 // 10**6
    subsec_a = subsec >> 8
    subsec_b = subsec & 0xFF
    uuid_int = (timestamp_ms & 0xFFFFFFFFFFFF) << 80
    uuid_int |= subsec_a << 64
    uuid_int |= subsec_b << 54
    uuid_int |= secrets.randbits(54)
    return uuid_int


def dict_to_formatted_str(d: Mapping[AnyStr, Any]) -> str:
    """Converts a dictionary into a formatted, JSON-like string.

    Align keys and values, including the alignment of sub-dicts.

    e.g.,
        >>> measurements = {'Hat': 'Large', 'Pants': {'Waist': 40, 'Inseam': 38}}
        >>> dict_to_formatted_str(measurements)
        Hat   : Large
        Pants : Waist  : 40
                Inseam : 38

    Args:
        d (dict): A dictionary with string-like keys.
    """
    # Initialize the output string
    output = ""
    # Get the system line separator
    nl = os.linesep
    # Determine the longest key
    key_len = 0
    for k in d.keys():
        key_len = max(len(k), key_len)
    # Add a space to the key length for clearance
    key_len += 1
    # Determine the tab space for sub-dictionaries
    tab = " " * (key_len + 2)
    # Iterate over the dictionary
    for k, v in d.items():
        # If the value is a dictionary itself...
        if isinstance(v, dict):
            # Then we write the key
            output += f"{k:<{key_len}}: "
            # Then call this function recursively on the value, and indent it
            # appropriately.
            output += tab.join(dict_to_formatted_str(v).splitlines(True))
        # Otherwise we just print the key and its value.
        elif isinstance(v, (int, float, complex)):
            output += f"{k:<{key_len}}: {v:,}" + nl
        else:
            output += f"{k:<{key_len}}: {v}" + nl
    return output


def check_current_folder(path, code):
    split_path = os.path.split(path)
    if split_path[-1] == "":
        split_path = os.path.split(split_path[0])
    return split_path[-1] == code
