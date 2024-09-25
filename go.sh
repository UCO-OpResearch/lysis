#!/bin/bash

exp_code="2024-04-16-1900"
in_code="_PLG2_tPA01_TB-xi.dat"
out_code="_TB-xi__1_582_867"
fort_executable="bin/macro_diffuse_into_and_along__external"

homedir="/home/bpaynter/git/UCO-OpResearch/lysis"
workdir="${homedir}"

cd "${homedir}/src/python"
python -u fortran_go.py     --in_code ${in_code}     --out_code ${out_code}.dat     -n 0     --cwd ${workdir}     ${workdir}/${fort_executable}     ${run_code}