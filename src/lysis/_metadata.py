"""Single source of truth for package metadata.

``__version__`` is resolved at install time by ``setuptools_scm`` (derived
from the nearest git tag) and read back here from the installed
distribution metadata.  ``__author__`` and ``__copyright__`` are the one
canonical place these strings are defined; everything else (the package
namespace, the CLI ``--version``, the Sphinx config) re-exports from here.
"""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("lysis")
except PackageNotFoundError:
    # Running from a raw source tree where the package was never installed
    # (so no dist metadata exists).  setuptools_scm only populates the
    # version at build/install time, so fall back to a sentinel.
    __version__ = "0.0.0+unknown"

__author__ = (
    "Brittany Bannish & Brad Paynter "
    "(see the AUTHORS file for the full list of contributors)"
)

__copyright__ = (
    "Copyright 2022-2026, Brittany Bannish <bbannish@uco.edu> "
    "& Brad Paynter <bpaynter@uco.edu>"
)
