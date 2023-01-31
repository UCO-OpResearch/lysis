#!/bin/bash -l
#SBATCH --job-name lysis2023-01-30-2002
#SBATCH --out /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2002/job.slurm-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

module purge
module load oneapi/compiler
/home/bpaynter/git/UCO-OpResearch/lysis/bin/macro_Q2_diffuse_into_and_along --expCode 2023-01-30-2002 --inFileCode _PLG2_tPA01_Q2.dat --outFileCode _PLG2_tPA01_into_and_along_Q2.dat > /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2002/macro_PLG2_tPA01_into_and_along_Q2.txt
  