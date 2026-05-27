# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
HOW TO MAINTAIN THIS FILE

- Every pull request that changes observable behaviour adds a bullet under
  "## [Unreleased]" in the appropriate group:
      Added       - new features
      Changed     - changes to existing behaviour
      Deprecated  - features slated for removal
      Removed      - features removed in this release
      Fixed        - bug fixes
      Security     - vulnerability fixes
- Write entries for users, not from commit messages. Reference the issue or
  PR number, e.g. "(#38)".
- At release time:
      1. Rename "## [Unreleased]" to "## [X.Y.Z] - YYYY-MM-DD".
      2. Add a fresh, empty "## [Unreleased]" section above it.
      3. Update the comparison links at the bottom of the file.
      4. Commit, then `git tag -a vX.Y.Z` (setuptools_scm derives the package
         version from the tag, so the tag is the source of truth).
- Semantic Versioning for the version number:
      MAJOR - backwards-incompatible change to the Python API, the CLI, or the
              on-disk data format.
      MINOR - new, backwards-compatible functionality.
      PATCH - backwards-compatible bug fix.
- Omit groups that have no entries. Keep the newest release at the top.
-->

## [Unreleased]

### Added

- This changelog, following the Keep a Changelog format. Releases v0.1.0–v0.3.0
  are summarised below from their git tags; finer-grained history before this
  file existed lives in the git log and the GitHub Releases.
- Command-line interface: the `lysis` command, with subcommands
  `init-experiment`, `init-macroscale`, `run-micro`, `run-macro`, `convert`,
  `validate`, `parameters`, `micro-stats`, `macro-stats`, `deg-rate`,
  `deg-time`, `compare`, `diff`, and `rename`.
- `analysis` subpackage: degradation analysis (fraction, rates, markers, lysis
  front), run/folder comparison, and HDF5 dataset diffing.
- Build and run provenance: each Fortran binary embeds its build commit and now
  its nearest release tag, reported via `--version` (#38); a preflight check
  refuses to run a stale binary (overridable with `--allow-stale-binary` or
  `LYSIS_ALLOW_STALE_BINARY=1`); binary and execution provenance are stamped
  onto the per-scale HDF5 groups; and `src/lysis/` commit-match gates guard the
  CLI run commands.
- Historical Fortran builds: rebuild and run a binary from an arbitrary past
  commit via `--fortran-commit`, with synthesised provenance.
- Packaging and tooling: `uv`-managed dependencies, `pyproject.toml`, and a
  pinned Python (`.python-version`); the package version is derived from the
  nearest git tag via `setuptools_scm` (#38); the license is declared in the
  package metadata (PEP 639, `GPL-3.0-only`) alongside the canonical `LICENSE`
  file (#38); and a single source of truth for `__version__`, `__author__`,
  and `__copyright__` lives in `src/lysis/_metadata.py` (#38).
- Continuous integration: a GitHub Actions test workflow.
- Documentation: system-onboarding, microscale-run, experiment-init, and
  historical-Fortran usage guides, an architecture overview, per-subpackage API
  reference pages, and an `AUTHORS` file.

### Changed

- Restructured the package: moved `src/python/lysis/` to a standard
  `src/lysis/` layout and split the monolithic `util/` submodule into focused
  subpackages (`config`, `dataio`, `geometry`, `execution`, `tools`).
- HDF5 data layer: macroscale output is written in the v2.0.0 specification,
  format conversion routes automatically through multi-step paths,
  macroscale input is generated lazily, and `DataStore.import_collection()`
  populates empty datasets from external sources.
- The CLI `--version` and the Sphinx `version`/`release`/`copyright`/`author`
  values are now derived from the package metadata (#38).

### Removed

- Per-file `__version__`, `__copyright__`, and `__author__` constants, plus the
  unused `__credits__`/`__license__`/`__maintainer__`/`__email__`/`__status__`
  header dunders, across all package modules (#38).

## [0.3.0] - 2026-05-19

### Added

- Python wrapper for the Fortran microscale model.
- Project ontology and terminology documentation.

## [0.2.0] - 2026-02-03

### Added

- Protofibril packing paper code.

## [0.1.0] - 2024-01-12

### Added

- Initial tagged release: the 2023 papers code (microscale and macroscale
  Fortran fibrinolysis model).

[Unreleased]: https://github.com/UCO-OpResearch/lysis/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/UCO-OpResearch/lysis/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/UCO-OpResearch/lysis/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/UCO-OpResearch/lysis/releases/tag/v0.1.0
