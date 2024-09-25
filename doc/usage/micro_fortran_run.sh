#!/bin/bash
#SBATCH --job-name=micro_rates
#SBATCH --output=micro_rates_%j.out
#SBATCH --nodes=1

### Usage
# 1. Fill out the folders in Setup
# 2. Add arguments at bottom of this file
# 3. Run this job
# 4. Examine micro_rates_###.out and micro______.txt

### Load Modules
module purge
module load oneapi/compiler

### Setup
# Replace these values MAKE SURE THERE ARE NO SPACES
LYSIS_ROOT=
RUN_CODE=


### Run
# eg. python3 your_script
cd $LYSIS_ROOT
mkdir -p data/$RUN_CODE
make

# Put parameters here on their own lines between --outFileCode and > data
# The format should be --param_name param_value \
bin/micro_rates \
    --runCode $RUN_CODE \
    --outFileCode $RUN_CODE.dat \
    --nodes 5 \
    > data/$RUN_CODE/micro_$RUN_CODE.txt

