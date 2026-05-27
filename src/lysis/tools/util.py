"""Utility functions for handling run data and parameters.

This module provides utility functions for managing experimental run data,
including formatting dictionaries into readable string representations.

Example usage:
    >>> from lysis.tools.util import dict_to_formatted_str
    >>> measurements = {'Hat': 'Large', 'Pants': {'Waist': 40, 'Inseam': 38}}
    >>> print(dict_to_formatted_str(measurements))
    Hat   : Large
    Pants : Waist  : 40
            Inseam : 38
"""

import os

from typing import Any, AnyStr, Mapping


def dict_to_formatted_str(d: Mapping[AnyStr, Any]) -> str:
    """Converts a dictionary into a formatted, JSON-like string.

    Aligns keys and values, including the alignment of sub-dictionaries.

    Example:
        >>> measurements = {'Hat': 'Large', 'Pants': {'Waist': 40, 'Inseam': 38}}
        >>> dict_to_formatted_str(measurements)
        Hat   : Large
        Pants : Waist  : 40
                Inseam : 38

    :param d: A dictionary with string-like keys
    :type d: dict
    :return: A formatted string representation of the dictionary
    :rtype: str
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
