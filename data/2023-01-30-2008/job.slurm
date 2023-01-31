#!/bin/bash -l
#SBATCH --job-name lysis2023-01-30-2008
#SBATCH --out /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2008/job.slurm-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

module purge
module load oneapi/compiler
/home/bpaynter/git/UCO-OpResearch/lysis/bin/macro_Q2_diffuse_into_and_along_fixed --expCode 2023-01-30-2008 --inFileCode _PLG2_tPA01_Kd00020036_Q2.dat --outFileCode _PLG2_tPA01_Kd00020036_into_and_along_fixed_Q2.dat --frac_forced 0.5143 --avgwait 277.8 --tf 1800 --seed -854989241 > /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2008/macro_PLG2_tPA01_Kd00020036_into_and_along_fixed_Q2.txt
  