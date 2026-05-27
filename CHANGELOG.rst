*********
Changelog
*********

All notable changes to this project are documented in this file.

The format is based on `Keep a Changelog
<https://keepachangelog.com/en/1.1.0/>`_, and this project adheres to
`Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

.. How to maintain this file:

   - Every pull request that changes observable behaviour adds a bullet under
     the "Unreleased" section, in the appropriate group:
         Added       - new features
         Changed     - changes to existing behaviour
         Deprecated  - features slated for removal
         Removed     - features removed in this release
         Fixed       - bug fixes
         Security    - vulnerability fixes
   - Write entries for users, not from commit messages. Reference the issue or
     PR number, e.g. "(#38)".
   - At release time:
       1. Rename the "Unreleased" heading to "X.Y.Z - YYYY-MM-DD".
       2. Add a fresh, empty "Unreleased" section above it.
       3. Update the links in the "Comparisons" section at the bottom.
       4. Commit, then `git tag -a vX.Y.Z` (setuptools_scm derives the package
          version from the tag, so the tag is the source of truth).
   - Semantic Versioning for the version number:
       MAJOR - backwards-incompatible change to the Python API, the CLI, or the
               on-disk data format.
       MINOR - new, backwards-compatible functionality.
       PATCH - backwards-compatible bug fix.
   - Omit groups that have no entries. Keep the newest release at the top.


Unreleased
==========

Added
-----

- This changelog, following the Keep a Changelog format. Releases v0.1.0–v0.3.0
  are summarised below from their git tags; finer-grained history before this
  file existed lives in the git log and the GitHub Releases.
- Command-line interface: the ``lysis`` command, with subcommands
  ``init-experiment``, ``init-macroscale``, ``run-micro``, ``run-macro``,
  ``convert``, ``validate``, ``parameters``, ``micro-stats``, ``macro-stats``,
  ``deg-rate``, ``deg-time``, ``compare``, ``diff``, and ``rename``.
- ``analysis`` subpackage: degradation analysis (fraction, rates, markers,
  lysis front), run/folder comparison, and HDF5 dataset diffing.
- Build and run provenance: each Fortran binary embeds its build commit and now
  its nearest release tag, reported via ``--version`` (#38); a preflight check
  refuses to run a stale binary (overridable with ``--allow-stale-binary`` or
  ``LYSIS_ALLOW_STALE_BINARY=1``); binary and execution provenance are stamped
  onto the per-scale HDF5 groups; and ``src/lysis/`` commit-match gates guard
  the CLI run commands.
- Historical Fortran builds: rebuild and run a binary from an arbitrary past
  commit via ``--fortran-commit``, with synthesised provenance.
- Packaging and tooling: ``uv``-managed dependencies, ``pyproject.toml``, and a
  pinned Python (``.python-version``); the package version is derived from the
  nearest git tag via ``setuptools_scm`` (#38); the license is declared in the
  package metadata (PEP 639, ``GPL-3.0-only``) alongside the canonical
  ``LICENSE`` file (#38); and a single source of truth for ``__version__``,
  ``__author__``, and ``__copyright__`` lives in ``src/lysis/_metadata.py``
  (#38).
- Continuous integration: a GitHub Actions test workflow.
- Documentation: system-onboarding, microscale-run, experiment-init, and
  historical-Fortran usage guides, an architecture overview, per-subpackage API
  reference pages, and an ``AUTHORS`` file.

Changed
-------

- Restructured the package: moved ``src/python/lysis/`` to a standard
  ``src/lysis/`` layout and split the monolithic ``util/`` submodule into
  focused subpackages (``config``, ``dataio``, ``geometry``, ``execution``,
  ``tools``).
- HDF5 data layer: macroscale output is written in the v2.0.0 specification,
  format conversion routes automatically through multi-step paths, macroscale
  input is generated lazily, and ``DataStore.import_collection()`` populates
  empty datasets from external sources.
- The CLI ``--version`` and the Sphinx ``version``/``release``/``copyright``/
  ``author`` values are now derived from the package metadata (#38).

Removed
-------

- Per-file ``__version__``, ``__copyright__``, and ``__author__`` constants,
  plus the unused ``__credits__``/``__license__``/``__maintainer__``/
  ``__email__``/``__status__`` header dunders, across all package modules
  (#38).


0.3.0 - 2026-05-19
==================

Adds a Python wrapper for the Fortran **microscale** model and takes several
steps toward the packaged layout that lands in v1.0.0; the architecture is
still "wrapper plus notebooks". A development release between papers, with no
associated publication.

Added
-----

- ``FortranMicro`` (in ``src/python/lysis/util/codeutil.py``), a companion to
  ``FortranMacro`` that drives the compiled microscale binary from Python, plus
  new entry points ``src/python/fortran_exec.py`` and ``cp_exec.py``.
- Project ontology (``docs/usage/ontology.rst``) defining Run, Experiment,
  Scenario, Mechanism, and Simulation.
- The first formal HDF5 v2.0.0 on-disk data specification
  (``docs/usage/data_specification.rst``).
- Stand-alone Fortran usage guides (``fortran_microscale.rst``,
  ``fortran_macroscale.rst``) and helper scripts (``micro_fortran_run.sh``,
  ``macro_fortran_run.sh``, ``octave_setup.sh``, ``micro_to_macro.py``,
  ``exec.sh``).
- New ``lysis.util`` modules: ``dataconvert.py`` (including the
  ``convert_fiber_degrade_time`` 1.99.0 → 2.0.0 converter), ``dataspec.py``
  (the ``DataSpec`` class), and ``fileops.py``; ``run.py`` splits ``Run`` out
  from ``Parameters``.

Changed
-------

- HDF5 / ``DataStore``: reformatted layout for structured data
  (``tpa_bind_events``, ``fiber_degrade_time``); micro and macro parameters
  stored as HDF attributes; logs added to the HDF group; ``DataStore`` can read
  HDF5 and export "micro to macro" data.
- Naming convention cleanup, kept in sync across Python and Fortran:
  ``expCode`` → ``runCode``; ``runs``/``stats`` → ``simulations``;
  ``nodes_in_row`` → ``nodes_in_micro_row``; seed parameters → ``micro_seed`` /
  ``macro_seed``; ``dist`` → ``radius``.
- Fortran sources: large reformatting of the macroscale files plus the
  ``runCode`` / ``radius`` renames; ``micro_rates.f90`` gains a
  ``snap_proportion`` parameter and ``trim()`` cleanups.

Removed
-------

- Legacy macroscale variants moved to ``src/fortran/_Archive/``;
  ``macro_Q2.f90`` removed (now archived).
- The four ``2024-0[1-3]-* - F-Macro Multi-Array-Process`` notebooks (their
  ``-Run`` companions are retained, so the simulation inputs remain
  reproducible).

Fixed
-----

- Microscale ``Lat`` allocation bug (the Austin tip commits carried in by the
  merge).


0.2.0 - 2026-02-03
==================

State used for the 2025 protofibril-packing-density paper (Risman et al.,
2025, *Research and Practice in Thrombosis and Haemostasis*,
`doi:10.1016/j.rpth.2025.102708 <https://doi.org/10.1016/j.rpth.2025.102708>`_).
The architecture is unchanged from v0.1.0 — still a Python wrapper plus
notebooks — but the codebase grew substantially.

Added
-----

- A GPL v3 ``LICENSE``.
- ReadTheDocs / Sphinx documentation skeleton (``docs/source/``,
  ``docs/usage/``, ``.readthedocs.yaml``, ``docs/requirements.txt``).
- New scenario and utility notebooks, and new analysis in the existing
  notebooks: degradation heatmaps, linear interpolation of individual-fiber
  degradation, and ``frac_forced`` calculations.

Changed
-------

- Fortran macroscale
  (``macro_diffuse_into_and_along__{external,internal}.f90``): precision fixes
  for ``Fort_Macro`` output and adoption of the new long ``f_deg``
  fiber-degrade-time data format.
- ``micro_rates.f90``: rate tables unchanged; indentation reflow, the
  ``expCode`` → ``runCode`` rename, and documentation of the Q0–Q3 fiber
  parameters.

Removed
-------

- Bulk simulation data removed from the git repository (now hosted externally;
  the dataset for this release is archived on Zenodo,
  `doi:10.5281/zenodo.15151618 <https://doi.org/10.5281/zenodo.15151618>`_).


0.1.0 - 2024-01-12
==================

First tagged release: an initial Python wrapper around the Fortran
**macroscale** model, with run orchestration and post-run analysis driven from
Jupyter notebooks rather than a packaged library. Mirrors the legacy
``v2023.final`` tag. Captures the code state used for the 2023–24 papers —
Risman et al. (2024, *Scientific Reports*,
`doi:10.1038/s41598-024-52844-4 <https://doi.org/10.1038/s41598-024-52844-4>`_)
and Bannish et al. (2024, *Biophysical Journal*,
`doi:10.1016/j.bpj.2024.02.002 <https://doi.org/10.1016/j.bpj.2024.02.002>`_).

Added
-----

- Python wrapper driving the Fortran macroscale simulation, with execution
  orchestration and post-run analysis performed in Jupyter notebooks.


Comparisons
===========

- `Unreleased
  <https://github.com/UCO-OpResearch/lysis/compare/v0.3.0...HEAD>`_
- `0.3.0 <https://github.com/UCO-OpResearch/lysis/compare/v0.2.0...v0.3.0>`_
- `0.2.0 <https://github.com/UCO-OpResearch/lysis/compare/v0.1.0...v0.2.0>`_
- `0.1.0 <https://github.com/UCO-OpResearch/lysis/releases/tag/v0.1.0>`_
