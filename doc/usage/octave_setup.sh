#!/bin/bash

module load Anaconda3
conda init
source ~/.bashrc

conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y --name lysis python=3.11
conda activate lysis
conda install -y black jupyterlab jupyter-black matplotlib pandas 
pip install octave-kernel