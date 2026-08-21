#!/usr/bin/env bash
#
# ood_jupyter_lab.sh -- prepare the uv environment for a bare `jupyter-lab` launch
#
# Intended for Open OnDemand (or any launcher) that SOURCES an environment
# setup script and then invokes `jupyter-lab` itself with its own flags (the
# websocket/proxy configuration OOD needs). This script therefore does NOT
# start JupyterLab: it only loads the required LMod modules and provisions
# and activates the project's uv-managed virtual environment, so that the
# subsequent bare `jupyter-lab` call resolves to .venv/bin/jupyter-lab and
# has Octave (for the octave-kernel) and the Intel compilers on PATH.
#
# It MUST be sourced, not executed -- the PATH/VIRTUAL_ENV changes have to
# persist into the shell that later runs `jupyter-lab`:
#
#   source /path/to/scripts/ood_jupyter_lab.sh
#
# In an OOD interactive-app template, point the environment setup / "custom
# script" field at this file and let OOD append the jupyter-lab command.
#
# Because it is sourced, this script deliberately avoids `set -e`/`exit`,
# which would terminate the caller's shell; it uses `return` on failure.
#
# NOTE (0.1.x line): the notebook stack (jupyterlab, octave-kernel,
# jupyter-black) is part of the package's MAIN dependencies on this branch,
# so a bare `uv sync` provisions everything -- there is no `notebooks` extra.

# --- Guard: this must be sourced, or the env changes are useless. -----------
if ! (return 0 2>/dev/null); then
    echo "ood_jupyter_lab.sh: must be SOURCED, not executed." >&2
    echo "  Use:  source ${BASH_SOURCE[0]:-$0}" >&2
    exit 1
fi

# --- Resolve the repository root from this script's own location. -----------
# ${BASH_SOURCE[0]} is the script path even when sourced.
_lysis_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Load required LMod modules. --------------------------------------------
# Octave backs the octave-kernel; the Intel compilers are needed for the
# Fortran binaries; FFmpeg backs matplotlib's FFMpegWriter, which the
# analysis notebooks use to render the simulation animations. Done before
# venv activation so .venv/bin stays at the front of PATH for the bare
# `jupyter-lab` call.
if ! command -v module >/dev/null 2>&1; then
    echo "ood_jupyter_lab.sh: 'module' (LMod) not available; cannot load modules." >&2
    unset _lysis_repo_root
    return 1
fi

module purge
module load Octave/10.1.0-foss-2023a
module load intel-compilers/2023.1.0
module load FFmpeg/6.0-GCCcore-12.3.0

# --- Provision the uv env (notebook stack is in the main dependencies). ------
# Run in a subshell so the caller's working directory is left unchanged.
if ! ( cd "$_lysis_repo_root" && uv sync ); then
    echo "ood_jupyter_lab.sh: 'uv sync' failed." >&2
    unset _lysis_repo_root
    return 1
fi

# --- Activate the venv so a bare `jupyter-lab` is on PATH. -------------------
_lysis_activate="$_lysis_repo_root/.venv/bin/activate"
if [[ ! -f "$_lysis_activate" ]]; then
    echo "ood_jupyter_lab.sh: venv activate script not found at $_lysis_activate" >&2
    unset _lysis_repo_root _lysis_activate
    return 1
fi

# shellcheck disable=SC1090
source "$_lysis_activate"

echo "ood_jupyter_lab.sh: activated $VIRTUAL_ENV"
echo "ood_jupyter_lab.sh: jupyter-lab -> $(command -v jupyter-lab)"

unset _lysis_repo_root _lysis_activate
