# Array Job Parallelism for the Microscale Simulation

## Motivation

The Fortran microscale binary (`micro_rates`) runs `micro_simulations` independent
Monte Carlo simulations in series within a single process.  Because each simulation
is completely independent (no shared state between iterations), they are a natural
candidate for data-parallel execution: split the simulation count across N jobs, run
them concurrently, and concatenate the results.

The goal of this work is to expose that parallelism through:

1. A `--array N` flag on the `lysis run-micro` CLI command.
2. A local (serial) reference implementation (`FortranMicro.run_array_full`) for
   testing without a Slurm cluster.
3. A Slurm job-array path (`#SBATCH --array=0-N-1`) for production use.

---

## What Was Implemented

### Branch: `hdf5`, commits `bff1a65` through `eabbb2c`

All work described below is on the `hdf5` branch.  The starting point (the last
commit *before* array work began) is `fa1ec9f`.

#### `src/lysis/execution/codeutil.py` — `FortranMicro`

* Added `n_array_jobs: int = None` field.
* **Array mode** (`index` and `n_array_jobs` both set):
  * Simulation distribution: job `i` runs
    `micro_simulations // n_array_jobs` simulations, with the first
    `micro_simulations % n_array_jobs` jobs each receiving one extra.
  * Seed derivation:
    `numpy.random.SeedSequence(base_seed).generate_state(n_array_jobs)[index]`
    (statistically independent seeds, NOT bit-exact vs. single-process run).
  * Output file code suffix: `__{index:02}` appended automatically.
* Added `import_array_results(staging_dir, hdf5_path, n_jobs, ...)` static method:
  * Locates per-chunk binary files (e.g. `lysis__00.dat`, `lysis__01.dat`, …).
  * Concatenates them into a `merged/` subdirectory inside `staging_dir`.
  * Reads total `micro_simulations` from the HDF5 `micro_data` group attributes
    (not from per-chunk `params.json`, which records only the chunk count).
  * Delegates to `import_results` (with `keep_tmpdir=True`) then manages the
    `merged/` lifecycle itself.
* Added `run_array_full(hdf5_path, n_jobs, *, keep_tmpdir=False)` method:
  * Runs N child `FortranMicro` instances sequentially in a shared `tmpdir`
    placed alongside the HDF5 file (not in `/tmp`).
  * Calls `import_array_results` after all children complete.
  * Cleanup policy: remove `tmpdir` on success unless `keep_tmpdir=True`;
    always preserve on failure.
* Updated `from_hdf5` to accept and forward `n_array_jobs`.

#### `src/lysis/tools/slurm.py`

* Added `generate_micro_array_child_script(...)` — generates a single Bash script
  with `#SBATCH --array=0-N-1`.  Each task uses `$SLURM_ARRAY_TASK_ID` as its
  `index` argument to `FortranMicro.from_hdf5`.  Supports both single-tier
  (default) and two-tier (`fast_tmp_root`) storage layouts.
* Added `_ARRAY_MASTER_PY_TEMPLATE` and `_generate_array_master_py(...)` — the
  master Python script submits the child array job via `gs.sbatch`, then polls
  `gs.squeue.read()` matching by base job-ID prefix until all tasks finish, then
  calls `FortranMicro.import_array_results`.
* Updated `submit_micro_slurm_job` to accept `n_array_jobs: Optional[int] = None`:
  * When `None` (default): existing single-child-script path unchanged.
  * When set (≥ 1): generates `child_array.sh` instead of `child_000.sh`, and
    uses the array master template.
  * Raises `ValueError` for `n_array_jobs < 1`.

#### `src/lysis/cli/run_micro.py`

* Added `--array N` option (`click.IntRange(min=1)`).
* With `--slurm`: forwards `n_array_jobs` to `submit_micro_slurm_job`; output
  message includes `(N array tasks)` alongside the job ID.
* Without `--slurm`: calls `fm.run_array_full(hdf5_path, n_jobs=N)`.

#### `docs/source/usage/run_microscale.rst`

* Added `--array N` option description with SeedSequence note.
* Added CLI examples: Slurm array, combined with partition/staging, local serial.
* Added `run_array_full` Python API section.
* Extended Typical Workflow to show array variant of the cluster loop.
* Added two troubleshooting entries: failed array tasks, wrong result count.

#### Tests

* `tests/execution/test_codeutil.py`:
  * `TestFortranMicroArrayExecCommand` — simulation distribution, seeds,
    remainder handling, file code suffix, reproducibility.
  * `TestFortranMicroFromHdf5Array` — `n_array_jobs` forwarding.
  * `TestFortranMicroImportArrayResults` — binary concatenation, merged
    `params.json`, cleanup.
  * `TestFortranMicroRunArrayFull` — N `exec_in_workdir` calls, correct child
    indices, `import_array_results` call, tmpdir location/cleanup.
  * `TestFortranArrayExecution` (`@pytest.mark.fortran_binary`) — runs
    `run_array_full` with 10 jobs × 100 simulations; checks dataset presence,
    total count = 1000, reproducibility, and statistical consistency vs.
    single-process run (within 15% tolerance).
* `tests/tools/test_slurm.py`:
  * `TestGenerateMicroArrayChildScript` — `--array` directive, `SLURM_ARRAY_TASK_ID`,
    `n_array_jobs` value, partition/`fast_tmp_root` forwarding.
  * `TestSubmitMicroSlurmJobArray` — end-to-end staging, correct script selection,
    master.py content, error on `n_array_jobs=0`.
* `tests/cli/test_run_micro.py`:
  * `TestRunMicroArrayLocal` — `run_array_full` dispatch, `n_jobs` forwarding,
    `keep_tmpdir`, `run_full` NOT called, `--array 0` rejected.
  * `TestRunMicroArraySlurm` — `submit_micro_slurm_job` called, `n_array_jobs`
    forwarded, job ID and task count in output, `n_array_jobs=None` when no `--array`.

---

## Known Limitation: Reproducibility

### The Problem

The current implementation uses `numpy.random.SeedSequence` to derive independent
per-chunk seeds from the base `micro_seed`.  Each chunk therefore starts from a
**different** RNG state than a single-process run.  The aggregate statistics (mean
lysis time, degradation fraction, etc.) are statistically equivalent, but the
results are **not bit-for-bit identical** to a single-process run of the same total
simulation count.

Consequently, the integration test (`TestFortranArrayExecution`) uses a 15%
statistical tolerance rather than an exact comparison, which is fragile and provides
weak validation.

### The Desired Behaviour

Running `N` array jobs of `micro_simulations / N` simulations should produce output
that is **bit-for-bit identical** to a single `micro_simulations` run with the same
base seed.

### The Required Change

The Fortran KISS32 RNG is a sequential stream seeded once at program start.  To
reproduce it across array jobs:

1. Job 0 starts from `base_seed` (normal initialization).
2. Job 0 saves its final KISS32 state (4 uint32 values: `c`, `jsr`, `x`, `y`) to a
   checkpoint file after completing its chunk.
3. Job 1 restores from that checkpoint file and continues the stream.
4. Repeat for all N jobs.

This requires the Fortran binary to support two new command-line flags:

* `--rngCheckpointIn FILE` — read 4 integers from `FILE` and call `set_kiss32`
  before the simulation loop (overrides seed-based initialization).
* `--rngCheckpointOut FILE` — call `get_kiss32` and write 4 integers to `FILE`
  after the simulation loop.

`kiss.c` already exposes `get_kiss32_` and `set_kiss32_` for exactly this purpose.

### Partial Implementation

A partial implementation of the Fortran changes was started (committed to this
branch as part of the pre-revert state):

```fortran
! New variables (after `integer :: seed = 0`):
character(80) :: rngCheckpointIn = ''
character(80) :: rngCheckpointOut = ''

! New command-line cases (in the select-case block):
case ('rngCheckpointIn')
    rngCheckpointIn = trim(param_value)
case ('rngCheckpointOut')
    rngCheckpointOut = trim(param_value)

! After call get_kiss32(stater):
if (len_trim(rngCheckpointIn) > 0) then
    open(99, file=trim(rngCheckpointIn), ...)
    read(99, *) stater(1), stater(2), stater(3), stater(4)
    close(99)
    call set_kiss32(stater)
    call get_kiss32(stater)
end if

! After end do  !enddo for stats loop:
if (len_trim(rngCheckpointOut) > 0) then
    call get_kiss32(stater)
    open(99, file=trim(rngCheckpointOut), ...)
    write(99, *) stater(1), stater(2), stater(3), stater(4)
    close(99)
end if
```

### Required Python-Side Changes

Once the Fortran binary supports checkpoints, `FortranMicro` needs:

1. New fields: `rng_checkpoint_in: str = None`, `rng_checkpoint_out: str = None`.
2. In `exec_command`: when `rng_checkpoint_in` is set, zero out `params["micro_seed"]`
   (so `--seed` is omitted) and append `--rngCheckpointIn FILE` instead.
   Append `--rngCheckpointOut FILE` when `rng_checkpoint_out` is set.
3. In `run_array_full`: chain checkpoints through the sequential jobs:
   * Job 0: `rng_checkpoint_in=None`, `rng_checkpoint_out=tmpdir/rng_00.txt`
   * Job i > 0: `rng_checkpoint_in=tmpdir/rng_{i-1:02d}.txt`,
     `rng_checkpoint_out=tmpdir/rng_{i:02d}.txt`
   * Remove the SeedSequence-based seed derivation from `run_array_full`'s child
     construction.
4. **Slurm path** (truly parallel): for Slurm the sequential checkpoint chain is not
   practical; the Slurm path should keep the SeedSequence approach OR adopt a
   pre-computation strategy (run N sequential seed-finding steps on the login node
   before submitting the array, saving N checkpoint files to the staging directory).

### Integration Test Change

With exact reproducibility, `TestFortranArrayExecution` can be simplified:
replace the three statistical consistency tests with a single `test_matches_single_process_exactly`
fixture that runs both paths and asserts bit-exact equality on all datasets.

---

## Files Changed

| File | Change |
|---|---|
| `src/lysis/execution/codeutil.py` | Array fields, `import_array_results`, `run_array_full` |
| `src/lysis/tools/slurm.py` | Array child script, array master, `n_array_jobs` param |
| `src/lysis/cli/run_micro.py` | `--array N` option |
| `src/fortran/micro_rates.f90` | Partial: `--rngCheckpointIn/Out` (not yet wired into Python) |
| `docs/source/usage/run_microscale.rst` | Array option docs, examples, troubleshooting |
| `tests/execution/test_codeutil.py` | Array unit tests + integration test |
| `tests/tools/test_slurm.py` | Slurm array unit tests |
| `tests/cli/test_run_micro.py` | CLI array unit tests |
