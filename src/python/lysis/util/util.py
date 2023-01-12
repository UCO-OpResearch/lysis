import os
from typing import Any, AnyStr, Mapping

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


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
    output = ''
    # Get the system line separator
    nl = os.linesep
    # Determine the longest key
    key_len = 0
    for k in d.keys():
        key_len = max(len(k), key_len)
    # Add a space to the key length for clearance
    key_len += 1
    # Determine the tab space for sub-dictionaries
    tab = ' ' * (key_len + 2)
    # Iterate over the dictionary
    for k, v in d.items():
        # If the value is a dictionary itself...
        if isinstance(v, dict):
            # Then we write the key
            output += f'{k:<{key_len}}: '
            # Then call this function recursively on the value, and indent it appropriately.
            output += tab.join(dict_to_formatted_str(v).splitlines(True))
        # Otherwise we just print the key and its value.
        elif isinstance(v, (int, float, complex)):
            output += f'{k:<{key_len}}: {v:,}' + nl
        else:
            output += f'{k:<{key_len}}: {v}' + nl
    return output
