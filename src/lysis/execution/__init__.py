from .run import *


def __getattr__(name):
    """Lazy import codeutil to avoid circular import with geometry.edge_grid."""
    from . import codeutil

    try:
        val = getattr(codeutil, name)
    except AttributeError:
        raise AttributeError(
            f"module {__name__!r} has no attribute {name!r}"
        )
    globals()[name] = val
    return val
