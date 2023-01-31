#!/bin/bash -l
#SBATCH --job-name lysis2023-01-30-2009
#SBATCH --out /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2009/job.slurm-%j.out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

module purge
module load oneapi/compiler
/home/bpaynter/git/UCO-OpResearch/lysis/bin/macro_Q2_diffuse_into --expCode 2023-01-30-2009 --inFileCode _PLG2_tPA01_Kd00020036_Q2.dat --outFileCode _PLG2_tPA01_Kd00020036_into_Q2.dat --frac_forced 0.5143 --avgwait 277.8 --seed -850336215 > /home/bpaynter/git/UCO-OpResearch/lysis/data/2023-01-30-2009/macro_PLG2_tPA01_Kd00020036_into_Q2.txt
  