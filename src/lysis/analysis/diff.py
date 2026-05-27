"""Element-wise data-table comparison between two simulation Runs.

Reads every non-log HDF5 dataset from each Run and compares the arrays
element-wise with an N-ULP floating-point tolerance (see
:func:`lysis.analysis.compare._values_match`).  Which scales
(``microscale_out``, ``macroscale_out``) are compared is decided by
:func:`available_scales` — the function auto-detects the presence of
each scale in the two Runs and compares the intersection.  The
asymmetric differences are reported back to the caller so the CLI can
warn when one Run lacks a scale the other has.

Typical workflow::

    from lysis.analysis.diff import diff_runs

    result = diff_runs(run1, run2)
    for label, entry in result["data_diff"].items():
        print(label, entry["status"])
    for scale in result["scales_skipped_in_run2"]:
        print(f"Warning: {scale} missing from run2")
"""


from typing import TYPE_CHECKING, Dict, Iterable, Set

import numpy as np

from lysis.analysis.compare import _compare_arrays

if TYPE_CHECKING:
    from lysis.config.run import Run

__all__ = [
    "available_scales",
    "diff_runs",
    "list_tables",
]

#: Dataset names (log tables) that must never be compared.
_LOG_DATASETS = frozenset({"micro_log", "macro_log"})

#: Scales that :func:`diff_runs` knows how to compare.
_KNOWN_SCALES = ("microscale_out", "macroscale_out")


def available_scales(run: "Run") -> Set[str]:
    """Return the subset of :data:`_KNOWN_SCALES` present in ``run``.

    Scales are detected via :attr:`~lysis.dataio.datastore.DataStore.collections`
    on the Run's open data store; each key in the returned set names a
    data collection the Run actually carries on disk.

    :param run: Run with data open.
    :type run: Run
    :return: Subset of ``{"microscale_out", "macroscale_out"}``.
    :rtype: set[str]
    """
    collections = run.data.collections
    return {scale for scale in _KNOWN_SCALES if scale in collections}


def _extract_data_arrays(
    run: "Run", scales: Iterable[str]
) -> Dict[str, np.ndarray]:
    """Extract every on-disk (non-log) data table from a Run's DataStore.

    Dataset contents are read into numpy arrays keyed by a display label
    of the form ``microscale_out/<name>`` or
    ``macroscale_out[<sim>]/<name>``.  Log tables (:data:`_LOG_DATASETS`)
    are skipped.  Derived datasets with ``data_location=None`` are
    already excluded by the :attr:`datasets` property on
    :class:`DataCollection` / :class:`SimulationView`.

    :param run: Run with data open.
    :type run: Run
    :param scales: Iterable subset of :data:`_KNOWN_SCALES` naming the
        collections to extract from.  Scales not present in the Run are
        silently skipped.
    :type scales: collections.abc.Iterable[str]
    :return: Mapping from label to numpy array copy.
    :rtype: dict[str, numpy.ndarray]
    """
    tables: Dict[str, np.ndarray] = {}
    data = run.data
    wanted = set(scales)
    collections = data.collections

    if "microscale_out" in wanted and "microscale_out" in collections:
        micro = data.microscale_out
        for name in micro.datasets:
            if name in _LOG_DATASETS:
                continue
            tables[f"microscale_out/{name}"] = getattr(micro, name)[:]

    if "macroscale_out" in wanted and "macroscale_out" in collections:
        macro = data.macroscale_out
        n_sims = run.macro_params.macro_simulations
        for sim in range(n_sims):
            view = macro[sim]
            for name in view.datasets:
                if name in _LOG_DATASETS:
                    continue
                tables[f"macroscale_out[{sim:02}]/{name}"] = getattr(view, name)[:]

    return tables


def list_tables(run: "Run") -> list:
    """Return a sorted list of every comparable dataset label in ``run``.

    Convenience wrapper used by the CLI to print an "Available tables:"
    list when ``--table`` does not match any dataset.

    :param run: Run with data open.
    :type run: Run
    :return: Sorted dataset labels in the same ``microscale_out/<name>``
        / ``macroscale_out[<sim>]/<name>`` form used by :func:`diff_runs`.
    :rtype: list[str]
    """
    return sorted(_extract_data_arrays(run, available_scales(run)))


def diff_runs(run1: "Run", run2: "Run") -> dict:
    """Compare on-disk data tables between two Runs element-wise.

    Auto-detects which scales each Run carries via
    :func:`available_scales` and compares the intersection.  If either
    Run lacks a scale the other has, the missing scale is listed in the
    returned dict so the caller can emit a warning.  Raises
    :class:`ValueError` if the two Runs have no scale in common (i.e.
    nothing to compare).

    For each dataset in the shared scales, reads the arrays from each
    Run and compares them element-wise.  Floating-point positions that
    differ by at most :data:`lysis.analysis.compare._MAX_ULPS` ULPs are
    treated as equivalent; every other dtype uses exact equality.
    Datasets that exist in one Run's scale but not the other's are
    reported with ``status="missing"``.

    :param run1: First Run with data open.
    :type run1: Run
    :param run2: Second Run with data open.
    :type run2: Run
    :return: Dict with four keys:

        - ``"data_diff"`` — ``{label: entry}`` mapping every dataset
          label to a result dict in the same format produced by
          :func:`lysis.analysis.compare._compare_arrays`.
        - ``"scales_compared"`` — sorted list of scales present in
          both Runs.
        - ``"scales_skipped_in_run1"`` — sorted list of scales present
          in ``run2`` but missing from ``run1``.
        - ``"scales_skipped_in_run2"`` — sorted list of scales present
          in ``run1`` but missing from ``run2``.
    :rtype: dict
    :raises ValueError: If neither Run has any of :data:`_KNOWN_SCALES`
        in common (no comparable data).
    """
    scales1 = available_scales(run1)
    scales2 = available_scales(run2)
    common = scales1 & scales2
    if not common:
        raise ValueError(
            "No comparable data collections found: neither run has any of "
            f"{sorted(_KNOWN_SCALES)} in common."
        )

    tables1 = _extract_data_arrays(run1, common)
    tables2 = _extract_data_arrays(run2, common)

    all_labels = sorted(set(tables1) | set(tables2))
    results: Dict[str, dict] = {}
    for label in all_labels:
        if label not in tables1:
            results[label] = {
                "status": "missing",
                "max_pct_diff": None,
                "location": None,
                "detail": "not in run1",
            }
        elif label not in tables2:
            results[label] = {
                "status": "missing",
                "max_pct_diff": None,
                "location": None,
                "detail": "not in run2",
            }
        else:
            results[label] = _compare_arrays(tables1[label], tables2[label])

    return {
        "data_diff": results,
        "scales_compared": sorted(common),
        "scales_skipped_in_run1": sorted(scales2 - scales1),
        "scales_skipped_in_run2": sorted(scales1 - scales2),
    }
