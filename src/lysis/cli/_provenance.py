"""Shared CLI helpers for the ``src/lysis/`` provenance gates.

Two pieces:

- :func:`allow_dirty_option` and :func:`allow_commit_mismatch_option` —
  Click-option decorators stacked onto each gated command.
- :func:`enforce_lysis_clean` and :func:`enforce_init_commit_match` —
  the actual gate logic, called from each command's body.  Both consume
  ``ctx.obj["lysis_dirty"]`` / ``ctx.obj["lysis_version"]`` populated
  by the top-level :func:`~lysis.cli.cli` group callback so they do not
  re-shell to git per command.

Centralising these here keeps the four CLI commands' integration to a
single line apiece.
"""

import warnings

import click

from ..config.constants import CONST
from ..tools.provenance import (
    allow_commit_mismatch_from_env,
    allow_dirty_from_env,
)


def allow_dirty_option(f):
    """Click decorator that adds ``--allow-dirty``.

    Downgrades the lysis-source dirty error to a warning.  Also honours
    the ``LYSIS_ALLOW_DIRTY`` environment variable.
    """
    return click.option(
        "--allow-dirty",
        "allow_dirty",
        is_flag=True,
        default=False,
        help=(
            "Proceed even if src/lysis/ has uncommitted changes.  The "
            "init/run stamp will still record the dirty state.  Also "
            f"honours {CONST.LYSIS_ALLOW_DIRTY_ENV}=1."
        ),
    )(f)


def allow_commit_mismatch_option(f):
    """Click decorator that adds ``--allow-commit-mismatch``.

    Used by ``run-micro`` / ``run-macro`` to bypass the check that the
    recorded init commit matches the currently-checked-out
    ``src/lysis/``.  Also honours ``LYSIS_ALLOW_COMMIT_MISMATCH=1``.
    """
    return click.option(
        "--allow-commit-mismatch",
        "allow_commit_mismatch",
        is_flag=True,
        default=False,
        help=(
            "Run even if the commit recorded by the matching init-* "
            "command differs from the currently-checked-out src/lysis/. "
            f"Also honours {CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV}=1."
        ),
    )(f)


def enforce_lysis_clean(ctx, allow_dirty):
    """Abort the current command if ``src/lysis/`` has uncommitted changes.

    Reads the cached dirty state from
    ``ctx.obj["lysis_dirty"]`` (populated by
    :func:`~lysis.cli.cli`).  Silent when clean or unknown; raises
    :class:`click.ClickException` when dirty and neither *allow_dirty*
    nor :func:`~lysis.tools.provenance.allow_dirty_from_env` permits
    proceeding.

    :param ctx: The Click context.
    :param allow_dirty: Truthy → override the error.
    :raises click.ClickException: When the gate fires.
    """
    if ctx.obj.get("lysis_dirty") != "dirty":
        return
    if allow_dirty or allow_dirty_from_env():
        return
    raise click.ClickException(
        "src/lysis/ has uncommitted changes; refusing to write provenance "
        "for an un-reproducible source tree.  Commit or stash your "
        "changes, or pass --allow-dirty (also honours "
        f"{CONST.LYSIS_ALLOW_DIRTY_ENV}=1) to proceed anyway."
    )


def enforce_init_commit_match(ctx, ds, scale, allow_commit_mismatch):
    """Verify the recorded init commit matches the current ``src/lysis/``.

    Compares ``ctx.obj["lysis_version"]`` against
    ``ds.read_init_provenance(scale)[init_version]``.  Behaviour:

    * Match → silent.
    * Mismatch + override allowed → :class:`UserWarning`, continue.
    * Mismatch + no override → :class:`click.ClickException`.
    * Legacy HDF5 file (no ``init_version`` attribute) →
      :class:`UserWarning` and continue.

    :param ctx: The Click context.
    :param ds: An open :class:`~lysis.dataio.datastore.DataStore`.
    :param scale: ``"micro"`` or ``"macro"``.
    :param allow_commit_mismatch: Truthy → downgrade mismatch error to
        a warning.
    :raises click.ClickException: On unoverridden mismatch.
    """
    recorded = ds.read_init_provenance(scale)
    current_version = ctx.obj.get("lysis_version", "unknown")

    if recorded is None:
        warnings.warn(
            f"HDF5 file has no '{CONST.INIT_VERSION_ATTR}' attribute "
            f"on its {scale}_params group — it pre-dates the lysis-"
            "source provenance feature; skipping commit-match check.",
            UserWarning,
            stacklevel=2,
        )
        return

    recorded_version = recorded.get(CONST.INIT_VERSION_ATTR, "unknown")
    if recorded_version == current_version:
        return

    override = allow_commit_mismatch or allow_commit_mismatch_from_env()
    msg = (
        f"src/lysis/ commit-match check failed for scale '{scale}': "
        f"init recorded {recorded_version!r} but the current src/lysis/ "
        f"resolves to {current_version!r}.  Check out the recorded "
        "commit, re-run the matching init-* command, or pass "
        "--allow-commit-mismatch (also honours "
        f"{CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV}=1)."
    )
    if not override:
        raise click.ClickException(msg)
    warnings.warn(msg, UserWarning, stacklevel=2)
