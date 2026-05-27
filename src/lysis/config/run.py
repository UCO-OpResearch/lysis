"""Container for run identity, parameters, and HDF5 data access.

This module provides the :class:`Run` class, which acts as the single entry
point for all information associated with a fibrinolysis simulation Run —
its identity (``run_code``), its on-disk layout, its micro- and macroscale
parameters, and a reference to its
:class:`~lysis.dataio.datastore.DataStore`.

Typical usage
-------------
Create a new Run and initialize parameters::

    from lysis.config.run import Run

    run = Run("/data/experiments")
    run.initialize_micro_param({"tpa_molecules": 500})
    run.initialize_macro_param({"rows": 64, "cols": 64})

Open the associated HDF5 DataStore for reading::

    with run.open_data() as ds:
        leaving_times = ds.microscale_out.tpa_leaving_time[:]
        pore_size = run.micro_params.pore_size
"""

import os
from datetime import datetime
from typing import Any, Mapping, Union

from ..dataio.datastore import DataStore
from ..tools.util import dict_to_formatted_str
from .parameters import MicroParameters, MacroParameters


__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class Run(object):
    """Houses all information about a given simulation Run.

    Holds the run's identity, on-disk paths, micro- and macroscale parameters,
    and a reference to the associated
    :class:`~lysis.dataio.datastore.DataStore`.  Parameters are initialised
    with :meth:`initialize_micro_param` and :meth:`initialize_macro_param`;
    the HDF5 data file is opened with :meth:`open_data`.

    **Instance attributes**

    .. attribute:: run_code
       :type: str

       Unique run identifier, used as the stem of the HDF5 filename
       (``{run_code}.h5``).  Auto-generated as ``YYYY-MM-DD-HHMM`` from
       the current date and time if not supplied.

    .. attribute:: os_path
       :type: str

       Absolute path to the directory containing the run's HDF5 file
       (i.e. ``data_root``).  The HDF5 file itself is at
       ``{os_path}/{run_code}.h5``.

    .. attribute:: micro_params
       :type: MicroParameters or None

       Microscale parameters.  ``None`` until :meth:`initialize_micro_param`
       is called.

    .. attribute:: macro_params
       :type: MacroParameters or None

       Macroscale parameters.  ``None`` until :meth:`initialize_macro_param`
       is called.  Requires ``micro_params`` to be set first.

    .. attribute:: data
       :type: DataStore or None

       Open :class:`~lysis.dataio.datastore.DataStore` for this run.
       ``None`` until :meth:`open_data` is called.

    :param data_root: Path to the directory that contains run HDF5 files.
    :type data_root: str, bytes, or os.PathLike
    :param run_code: Identifier for this Run.  Must be a date-time string in
        ``YYYY-MM-DD-HHMM`` format.  If omitted, one is generated from the
        current date and time.
    :type run_code: str, optional
    :raises RuntimeError: If ``data_root`` is not an existing directory.
    """

    def __init__(self, data_root: Union[str, bytes, os.PathLike], run_code: str = None):
        # Check if the data folder path is valid
        if not os.path.isdir(data_root):
            raise RuntimeError("Data folder not found.", data_root)
        self.os_path = str(data_root)
        # If no run code was given, create a new one from the current
        # date and time.
        if run_code is None:
            self.run_code = datetime.now().strftime("%Y-%m-%d-%H%M")
        else:
            self.run_code = run_code

        # TODO(bpaynter): Check if the parameters are already stored.
        #                 Don't allow parameters to be changed once stored.
        # Initialize the internal storage as empty
        # self.sequence: ExpComponent = ExpComponent.NONE
        self.micro_params = None
        self.macro_params = None
        self.data = None
        self._cache = {}

    def __str__(self) -> str:
        """Return a human-readable formatted string of the run's parameters.

        :return: Formatted representation of ``run_code``, ``micro_params``,
            and ``macro_params``.
        :rtype: str
        """
        # Convert internal storage to a dictionary
        values = self.to_dict()
        # Call the formatter and return
        return dict_to_formatted_str(values)

    def initialize_micro_param(self, params: Mapping[str, Any] = None) -> None:
        """Set the microscale parameters for this Run.

        Constructs a :class:`~lysis.config.parameters.MicroParameters` object
        using default values, overriding any keys supplied in ``params``.
        The result is stored as ``self.micro_params``.

        :param params: Parameter values that differ from the defaults,
            e.g. ``{"binding_rate": 10, "pore_size": 3}``.
        :type params: dict, optional
        """
        if params is not None:
            self.micro_params = MicroParameters(**params)
        else:
            self.micro_params = MicroParameters()

    def initialize_macro_param(self, params: dict[str, Any] = None) -> None:
        """Set the macroscale parameters for this Run.

        Constructs a :class:`~lysis.config.parameters.MacroParameters` object
        using default values, overriding any keys supplied in ``params``.
        The result is stored as ``self.macro_params``.

        :meth:`initialize_micro_param` must be called before this method,
        as macroscale parameters are derived from the microscale parameters.

        :param params: Parameter values that differ from the defaults,
            e.g. ``{"rows": 64, "cols": 64}``.
        :type params: dict, optional
        :raises RuntimeError: If ``micro_params`` has not been set.
        """
        # The macroscale model is dependent on the parameters and results of the
        # microscale model. If no microscale parameters are supplied, the macroscale
        # model cannot be initialized.
        if self.micro_params is None:
            raise RuntimeError("No Microscale parameters.")
        if params is not None:
            self.macro_params = MacroParameters(self.micro_params, **params)
        else:
            self.macro_params = MacroParameters(micro_params=self.micro_params)

    def load_params_from_hdf5(self) -> None:
        """Load micro- and macroscale parameters from the run's HDF5 file.

        Opens ``{run_code}.h5`` in read-only mode and populates
        ``self.micro_params`` from the stored microscale attributes.  If
        macroscale parameters are present in the file, ``self.macro_params``
        is populated as well.

        This is the inverse of creating a DataStore with
        :meth:`~lysis.dataio.datastore.DataStore.create` / storing params via
        ``initialize_micro_param``.  Use it to reconstruct a fully-parameterised
        :class:`Run` from an existing HDF5 file without knowing the original
        parameter values.

        :raises RuntimeError: If ``os_path`` does not contain ``{run_code}.h5``.
        :raises ValueError: If the HDF5 file's dataspec version is incompatible.
        """
        with DataStore(self.run_code, self.os_path, mode="r") as ds:
            self.micro_params = ds.micro_params
            if ds.macro_params is not None:
                self.macro_params = ds.macro_params

    def rename(self, new_run_code: str) -> str:
        """Rename this Run's HDF5 file and record the previous name.

        Delegates the on-disk rename and ``renamed_from`` history update
        to :meth:`lysis.dataio.datastore.DataStore.rename`, then updates
        this Run's ``run_code`` to reflect the new identifier.  Any
        currently open :class:`~lysis.dataio.datastore.DataStore` is
        closed first so the file handle does not point at the old path.

        :param new_run_code: The new run identifier (becomes the new HDF5
            filename stem).
        :type new_run_code: str
        :return: The previous ``run_code``.
        :rtype: str
        :raises ValueError: If ``new_run_code`` equals the current
            ``run_code``.
        :raises FileNotFoundError: If the current HDF5 file does not exist.
        :raises FileExistsError: If an HDF5 file already exists at the
            destination path.
        """
        old_run_code = self.run_code
        if self.data is not None:
            self.data.close()
            self.data = None
        DataStore.rename(old_run_code, new_run_code, self.os_path)
        self.run_code = new_run_code
        return old_run_code

    def open_data(self, mode: str = "r") -> DataStore:
        """Open the HDF5 DataStore for this run.

        Creates a :class:`~lysis.dataio.datastore.DataStore` backed by the HDF5
        file ``{run_code}.h5`` in the run's data directory. The DataStore is
        stored as ``self.data`` and also returned for convenience.

        :param mode: HDF5 file mode passed to :class:`DataStore`.
            ``"r"`` (default) for read-only, ``"a"`` for read/write.
        :type mode: str
        :return: The opened DataStore.
        :rtype: DataStore
        """
        self.data = DataStore(self.run_code, self.os_path, mode=mode)
        return self.data

    def to_dict(self) -> dict:
        """Return the run's parameters as a plain dictionary.

        Does not include system-specific information such as paths.  The
        returned dictionary has the following keys:

        - ``"run_code"`` — the run identifier string.
        - ``"micro_params"`` — base-dict of microscale parameters, or
          ``None`` if not yet initialised.
        - ``"macro_params"`` — base-dict of macroscale parameters, or
          ``None`` if not yet initialised.

        :return: Dictionary representation of the run's parameters.
        :rtype: dict
        """
        # Initialize a dictionary of the appropriate parameters
        output = {
            "run_code": self.run_code,
            "data_filenames": None,
            "micro_params": None,
            "macro_params": None,
        }
        # Get the data filenames from the DataStore
        # if self.data is not None:
        #     output["data_filenames"] = self.data.to_dict()
        # Convert the Microscale parameters to a dictionary
        if self.micro_params is not None:
            output["micro_params"] = self.micro_params.to_basedict()

        # Convert the Macroscale parameters to a dictionary
        if self.macro_params is not None:
            output["macro_params"] = self.macro_params.to_basedict()
        return output
