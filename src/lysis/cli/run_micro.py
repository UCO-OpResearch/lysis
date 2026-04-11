"""``lysis run-micro`` — execute the Fortran microscale simulation for a Run.

Reads micro-parameters from an existing HDF5 file, executes the compiled
Fortran microscale binary, and imports the results back into the same file.

Execution modes
---------------
**Local** (default)
    Creates a temporary working directory alongside the HDF5 file, runs the
    binary synchronously, imports results, and cleans up.

**Slurm** (``--slurm``)
    Submits a master Slurm job that orchestrates one or more child Fortran
    jobs, polls for completion, imports results, and optionally cleans up.
    Returns immediately after printing the master job ID.
"""

import click

from lysis.cli import cli


@cli.command(name="run-micro")
@click.argument("hdf5_path", metavar="HDF5_PATH", type=click.Path(exists=True))
@click.option(
    "--executable",
    required=True,
    type=click.Path(),
    help="Path to the compiled Fortran microscale binary.",
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
        "Defaults to the parent directory of HDF5_PATH.  "
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
        "When set, the Fortran binary writes here and results are moved "
        "to the staging directory on completion (two-tier storage).  "
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
    "--file-code",
    "file_code",
    default="",
    metavar="TEXT",
    show_default=True,
    help="Output file code suffix for the Fortran binary.",
)
@click.option(
    "--array",
    "n_array_jobs",
    default=None,
    type=click.IntRange(min=1),
    metavar="N",
    help=(
        "Split the simulation into N parallel array jobs.  "
        "With --slurm, submits a Slurm job array (``#SBATCH --array=0-N-1``); "
        "without --slurm, executes N sub-jobs serially in the same process "
        "(useful for testing).  Simulations are distributed evenly across "
        "jobs with independent seeds derived from the base micro_seed."
    ),
)
@click.pass_context
def run_micro(ctx, hdf5_path, executable, use_slurm, partition, staging_root,
              fast_tmp_root, keep_tmpdir, file_code, n_array_jobs):
    """Execute the Fortran microscale simulation for a Run.

    HDF5_PATH must point to an existing ``.h5`` file containing
    ``micro_params``.  The simulation output is imported back into that same
    file on completion.

    \b
    Examples:
        lysis run-micro data/run01.h5 --executable bin/micro.exe
        lysis run-micro data/run01.h5 --executable bin/micro.exe --keep-tmpdir
        lysis run-micro data/run01.h5 --executable bin/micro.exe --slurm
        lysis run-micro data/run01.h5 --executable bin/micro.exe --slurm \\
            --partition normal --staging-root /scratch/staging
        lysis run-micro data/run01.h5 --executable bin/micro.exe --slurm \\
            --fast-tmp-root /nvme/scratch
        lysis run-micro data/run01.h5 --executable bin/micro.exe --slurm \\
            --array 10
    """
    console = ctx.obj["console"]

    if use_slurm:
        from lysis.tools.slurm import submit_micro_slurm_job

        job_id = submit_micro_slurm_job(
            hdf5_path,
            executable,
            staging_root=staging_root,
            partition=partition,
            fast_tmp_root=fast_tmp_root,
            keep_tmpdir=keep_tmpdir,
            out_code=file_code,
            n_array_jobs=n_array_jobs,
        )
        if n_array_jobs is not None:
            console.print(
                f"Submitted master Slurm job [bold]{job_id}[/bold] "
                f"({n_array_jobs} array tasks)"
            )
        else:
            console.print(f"Submitted master Slurm job [bold]{job_id}[/bold]")
    else:
        from lysis.execution.codeutil import FortranMicro

        fm = FortranMicro.from_hdf5(hdf5_path, executable,
                                     out_file_code=file_code)
        if n_array_jobs is not None:
            if not ctx.obj.get("verbose", 0):
                with console.status(
                    f"Running {n_array_jobs}-job array for {fm.run.run_code}..."
                ):
                    fm.run_array_full(hdf5_path, n_jobs=n_array_jobs,
                                      keep_tmpdir=keep_tmpdir)
            else:
                console.print(
                    f"Running [bold]{n_array_jobs}[/bold]-job array for "
                    f"[bold]{fm.run.run_code}[/bold]"
                )
                fm.run_array_full(hdf5_path, n_jobs=n_array_jobs,
                                  keep_tmpdir=keep_tmpdir)
        else:
            if not ctx.obj.get("verbose", 0):
                with console.status(
                    f"Running microscale simulation for {fm.run.run_code}..."
                ):
                    fm.run_full(hdf5_path, keep_tmpdir=keep_tmpdir)
            else:
                console.print(
                    f"Running microscale simulation for [bold]{fm.run.run_code}[/bold]"
                )
                fm.run_full(hdf5_path, keep_tmpdir=keep_tmpdir)
        console.print(
            f"[green]Microscale results imported into[/green] {hdf5_path}"
        )
