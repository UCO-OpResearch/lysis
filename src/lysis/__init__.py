from ._metadata import __author__, __copyright__, __version__
from .molecule import *
from .np_macroscale import *

try:
    import cupy
    from .cp_macroscale import *
except (ImportError, TypeError):
    # ImportError: cupy not installed
    # TypeError: cupy is mocked (e.g. by Sphinx autodoc) and mock objects
    #   don't support type union syntax (cp.ndarray | None)
    pass
