"""``lysis run-macro`` — execute the Fortran macroscale simulation for a Run.

Reads micro- and macro-parameters from an existing HDF5 file, generates the
required macroscale input files, executes the compiled Fortran macroscale
binary for each simulation, and imports the results back into the same file.

Accepts either a single ``.h5`` file or a directory containing multiple HDF5
files (an experiment folder).  When a directory is given, all runs are
processed in order.

Execution modes
---------------
**Local** (default)
    Creates a temporary working directory alongside the HDF5 file, runs the
    binary synchronously, imports results, and cleans up.  In batch mode, runs
    are processed sequentially; the first failure stops the batch.

**Slurm** (``--slurm``)
    Submits a master Slurm job that pre-stages input files, orchestrates a
    Slurm array job (one task per simulation), imports results, and optionally
    cleans up.  Returns immediately after printing the master job ID(s).  In
    batch mode, all master jobs are submitted in parallel (without waiting for
    earlier jobs to finish).
"""

import contextlib
from pathlib import Path

import click

from lysis.cli import cli
from lysis.cli._provenance import (
    allow_commit_mismatch_option,
    allow_dirty_option,
    enforce_init_commit_match,
    enforce_lysis_clean,
)
from lysis.tools.slurm import DEFAULT_MODULES, parse_sbatch_tokens


@cli.command(name="run-macro")
@click.argument("target_path", metavar="PATH", type=click.Path(exists=True))
@click.option(
    "--executable",
    required=True,
    type=click.Path(),
    help="Path to the compiled Fortran macroscale binary.",
)
@click.option(
    "--slurm",
    "use_slurm",
    is_flag=True,
    default=False,
    help="Dispatch via Slurm.  Submits a master job and prints the job ID.",
)
@click.option(
    "--partition",
    default=None,
    metavar="NAME",
    help="Slurm partition (only meaningful with --slurm).",
)
@click.option(
    "--staging-root",
    "staging_root",
    default=None,
    type=click.Path(),
    metavar="PATH",
    help=(
        "Root directory for the shared staging temp dir.  "
        "Defaults to the parent directory of the HDF5 file.  "
        "Only meaningful with --slurm."
    ),
)
@click.option(
    "--fast-tmp-root",
    "fast_tmp_root",
    default=None,
    type=click.Path(),
    metavar="PATH",
    help=(
        "Root directory for fast node-local scratch storage.  "
        "When set, each Slurm array task writes here and moves results to "
        "the staging directory on completion (two-tier storage).  "
        "Only meaningful with --slurm."
    ),
)
@click.option(
    "--keep-tmpdir",
    "keep_tmpdir",
    is_flag=True,
    default=False,
    help="Always preserve the temporary output directory (useful for debugging).",
)
@click.option(
    "--out-file-code",
    "out_file_code",
    default="",
    metavar="TEXT",
    show_default=True,
    help=(
        "Output file code suffix for the Fortran binary.  "
        "Not supported when PATH is a directory."
    ),
)
@click.option(
    "--in-file-code",
    "in_file_code",
    default="",
    metavar="TEXT",
    show_default=True,
    help="Input file code suffix for macroscale_in files.",
)
@click.option(
    "--modules",
    "modules",
    default=DEFAULT_MODULES,
    show_default=True,
    metavar="MODULES",
    help=(
        "Space-separated list of LMod module specs to ``module load`` in "
        "each generated Slurm script (e.g. ``\"intel-compilers/2024 "
        "SciPy-bundle/2023.07\"``).  Include a Fortran compiler module so "
        "the simulation binary finds its runtime.  With --slurm, the same "
        "list also wraps the historical-build ``make`` invocation of "
        "--fortran-commit so the binary links against the matching "
        "toolchain.  Only meaningful with --slurm."
    ),
)
@click.option(
    "--fortran-commit",
    "fortran_commit",
    default=None,
    metavar="REF",
    help=(
        "Build the Fortran binary from this git ref (SHA, tag, or branch) "
        "before running, instead of using a pre-built binary at "
        "--executable's path.  When set, --executable is interpreted as the "
        "binary name (with or without a leading 'bin/') to pick from the "
        "historical build's bin/ directory; the binary↔src/fortran "
        "staleness check is bypassed because the mismatch is intentional.  "
        "Synthesised provenance (resolved SHA, ``iso_fortran_env`` compiler "
        "string, ``binary_source = 'historical:<sha>'``) is stamped into "
        "the HDF5 file in place of the binary's own --version output."
    ),
)
@click.option(
    "--sbatch",
    "sbatch_tokens",
    multiple=True,
    metavar="TOKEN",
    help=(
        "Override #SBATCH options on the generated master and array task "
        "scripts.  Repeatable.  Syntax: 'KEY=VALUE' sets/overrides, 'KEY' "
        "alone adds a flag-style option (e.g. 'hold'), and '^KEY' removes "
        "a default whose dict key matches KEY exactly (use this to drop "
        "the built-in 'exclusive=user' via --sbatch '^exclusive=user').  "
        "Splits on the first '='.  Only meaningful with --slurm."
    ),
)
@click.option(
    "--allow-stale-binary",
    "allow_stale_binary",
    is_flag=True,
    default=False,
    help=(
        "Run even if the Fortran binary's embedded build commit does not "
        "match the currently-checked-out src/fortran/.  A loud warning is "
        "written to stderr, prepended to each Fortran log file, and stamped "
        "onto the resulting HDF5 group.  Also honors the env var "
        "LYSIS_ALLOW_STALE_BINARY=1."
    ),
)
@allow_dirty_option
@allow_commit_mismatch_option
@click.pass_context
def run_macro(ctx, target_path, executable, use_slurm, partition, staging_root,
              fast_tmp_root, keep_tmpdir, out_file_code, in_file_code, modules,
              fortran_commit, sbatch_tokens, allow_stale_binary, allow_dirty,
              allow_commit_mismatch):
    """Execute the Fortran macroscale simulation for a Run or Experiment.

    PATH may be either:

    \b
      - A single ``.h5`` file containing both ``micro_params`` and
        ``macro_params``.
      - A directory (experiment folder).  If ``experiment.json`` is present,
        runs are processed in CSV order; otherwise all ``*.h5`` files in the
        directory are processed in alphabetical order.

    The simulation output is imported back into each HDF5 file on completion.

    \b
    Examples:
        lysis run-macro data/run01.h5 --executable bin/macro.exe
        lysis run-macro data/run01.h5 --executable bin/macro.exe --keep-tmpdir
        lysis run-macro data/run01.h5 --executable bin/macro.exe --slurm
        lysis run-macro data/run01.h5 --executable bin/macro.exe --slurm \\
            --partition normal --staging-root /scratch/staging
        lysis run-macro data/run01.h5 --executable bin/macro.exe --slurm \\
            --fast-tmp-root /nvme/scratch
        lysis run-macro data/my-experiment/ --executable bin/macro.exe
        lysis run-macro data/my-experiment/ --executable bin/macro.exe --slurm
    """
    console = ctx.obj["console"]

    # Gate: refuse to write provenance if src/lysis/ is dirty (unless overridden).
    enforce_lysis_clean(ctx, allow_dirty)

    target = Path(target_path)
    is_batch = target.is_dir()

    if is_batch:
        if out_file_code:
            raise click.UsageError(
                "--out-file-code is not supported when PATH is a directory."
            )
        from lysis.cli._batch import resolve_hdf5_paths
        hdf5_paths = resolve_hdf5_paths(target)
    else:
        hdf5_paths = [target]

    # Verify the init-macroscale commit recorded on each HDF5 file matches
    # the currently-checked-out src/lysis/ commit, before any work starts.
    from lysis.dataio.datastore import DataStore
    for hdf5_path in hdf5_paths:
        with DataStore(hdf5_path.stem, str(hdf5_path.parent), mode="r") as ds:
            enforce_init_commit_match(ctx, ds, "macro", allow_commit_mismatch)

    # --fortran-commit: build the historical binary once up front, share its
    # path across every run in this invocation, tear down at the end.  Only
    # load the modules around the build when Slurm jobs (which load the same
    # modules via the preamble) will be executing the binary.
    with contextlib.ExitStack() as stack:
        if fortran_commit is not None:
            from lysis.execution.historical_build import (
                HistoricalBuildError,
                build_historical_binary,
            )
            try:
                resolved_exe, historical_provenance = stack.enter_context(
                    build_historical_binary(
                        fortran_commit,
                        executable,
                        modules=modules if use_slurm else None,
                        keep_dir=keep_tmpdir,
                    )
                )
            except HistoricalBuildError as e:
                raise click.ClickException(str(e))
            executable = str(resolved_exe)
            sha = historical_provenance.get("binary_commit", "?")
            console.print(
                f"[yellow]Using Fortran built from {sha[:7]} "
                f"(--fortran-commit {fortran_commit}); "
                f"binary↔src/fortran mismatch is expected and not checked."
                f"[/yellow]"
            )
        else:
            historical_provenance = None

        if use_slurm:
            from lysis.tools.slurm import submit_macro_slurm_job

            try:
                sbatch_overrides = parse_sbatch_tokens(sbatch_tokens)
            except ValueError as e:
                raise click.BadParameter(str(e), param_hint="--sbatch")

            submitted = []
            for hdf5_path in hdf5_paths:
                try:
                    job_id = submit_macro_slurm_job(
                        hdf5_path,
                        executable,
                        staging_root=staging_root,
                        partition=partition,
                        fast_tmp_root=fast_tmp_root,
                        keep_tmpdir=keep_tmpdir,
                        in_code=in_file_code,
                        out_code=out_file_code,
                        modules=modules,
                        sbatch_overrides=sbatch_overrides,
                        historical_binary_attrs=historical_provenance,
                    )
                except ValueError as e:
                    raise click.ClickException(str(e))
                submitted.append((hdf5_path.stem, job_id))
                console.print(
                    f"Submitted master Slurm job [bold]{job_id}[/bold]"
                    f" for {hdf5_path.stem}"
                )

            if len(submitted) > 1:
                console.print(
                    f"[green]Submitted {len(submitted)} master Slurm jobs.[/green]"
                )
        else:
            from lysis.execution.fortran_macro import FortranMacro

            n = len(hdf5_paths)
            for i, hdf5_path in enumerate(hdf5_paths):
                prefix = f"[{i + 1}/{n}] " if n > 1 else ""

                try:
                    fm = FortranMacro.from_hdf5(
                        hdf5_path,
                        executable,
                        in_file_code=in_file_code,
                        out_file_code=out_file_code,
                        allow_stale_binary=allow_stale_binary,
                        skip_binary_verification=(
                            historical_provenance is not None
                        ),
                        historical_binary_attrs=historical_provenance,
                    )
                except ValueError as e:
                    raise click.ClickException(str(e))
                if not ctx.obj.get("verbose", 0):
                    with console.status(
                        f"{prefix}Running macroscale simulation for"
                        f" {fm.run.run_code}..."
                    ):
                        fm.run_full(keep_tmpdir=keep_tmpdir)
                else:
                    console.print(
                        f"{prefix}Running macroscale simulation for"
                        f" [bold]{fm.run.run_code}[/bold]"
                    )
                    fm.run_full(keep_tmpdir=keep_tmpdir)
                console.print(
                    f"[green]Macroscale results imported into[/green] {hdf5_path}"
                )
