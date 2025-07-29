#!/bin/bash

module load Miniforge
conda init
source ~/.bashrc

conda create -y --name lysis python=3.11
conda activate lysis
conda install -y black jupyterlab jupyter-black \
   matplotlib pandas numpy
pip install octave-kernel