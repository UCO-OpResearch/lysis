#!/bin/bash
#SBATCH --job-name=macro
#SBATCH --output=macro__%j.out
#SBATCH --nodes=1
#SBATCH --array=0
#SBATCH --exclusive=user 

### Usage
# 1. Fill out the folders in Setup
# 2. Add arguments at bottom of this file
# 3. Run this job
# 4. Examine macro__###.out and macro__XXXXXXX.txt

### Load Modules
module purge
module load intel-compilers/2023.1.0
module load SciPy-bundle/2023.07-gfbf-2023a

### Setup
SIM=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
# Replace these values MAKE SURE THERE ARE NO SPACES
LYSIS_ROOT=/home/bbannish/tPA_variant/lysis
MICRO_FILE_CODE_BIG=variant3p6
MICRO_FILE_CODE_SMALL=variant0p036
MACRO_RUN_CODE=macro_big3p6_small0p036

# Macroscale grid size
N=93
F=121

### Run
cd $LYSIS_ROOT
mkdir -p data/$MACRO_RUN_CODE/$SIM
make

# Produce the input files the Macroscale 
FRAC_FORCED_BIG=$(python3 doc/usage/micro_to_macro.py \
    --run_code $MACRO_RUN_CODE/$SIM \
    --in_code $MICRO_FILE_CODE_BIG \
    --out_code $MACRO_RUN_CODE \
    -N $N \
    -F $F)
    
FRAC_FORCED_SMALL=$(python3 doc/usage/micro_to_macro.py \
    --run_code $MACRO_RUN_CODE/$SIM \
    --in_code $MICRO_FILE_CODE_SMALL \
    --out_code $MACRO_RUN_CODE \
    -N $N \
    -F $F)

# Put parameters here on their own lines between --outFileCode and > data
# The format should be --param_name param_value \
bin/macro_diffuse_into_and_along_variant__external \
    --runCode $MACRO_RUN_CODE/$SIM \
    --inFileCode $MICRO_FILE_CODE \
    --outFileCode $MACRO_RUN_CODE \
    --N $N \
    --F $F \
    --frac_forced_big $FRAC_FORCED_BIG \
    --frac_forced_small $FRAC_FORCED_SMALL \
    --avgwait_big 2.78 \
    --avgwait_small 277.8 \
    > data/$MACRO_RUN_CODE/$SIM/macro__$MACRO_RUN_CODE.txt
