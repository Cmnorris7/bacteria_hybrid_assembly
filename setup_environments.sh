#!/bin/bash

# This script will clone Autocycler repo and setup conda environments

# Create qc_polish_hybrid conda environment
conda env create --file environment.yml -n qc_polish_hybrid 

# Clone Autocycler repo
cd ../ && git clone https://github.com/rrwick/Autocycler.git

# Create autocycler conda environment
cd ../Autocycler/scripts/
conda env create --file environment.yml --name autocycler 
conda activate autocycler
cp *.py *.sh "$CONDA_PREFIX"/bin  # copy assembly helper scripts into conda env
plassembler download -d "$CONDA_PREFIX"/plassembler_db  # download plassembler_db
conda create -f environment.yml ../