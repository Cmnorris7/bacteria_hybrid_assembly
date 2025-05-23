#!/bin/bash

# This script will clone Autocycler repo and setup conda environments

# Create qc_polish_hybrid conda environment
conda env create -y --file ./config/environment.yml -n qc_polish_hybrid 

# Clone Autocycler repo
cd ../ && git clone https://github.com/rrwick/Autocycler.git

# Create autocycler conda environment
cd Autocycler/scripts/
conda env create -y --file environment.yml --name autocycler

# Activate autocycler conda environment
CONDA_ROOT=$(conda info --base)
source "$CONDA_ROOT/etc/profile.d/conda.sh" && conda activate autocycler

# Copy helper scripts and download plassembler_db
cp *.py *.sh "$CONDA_PREFIX"/bin  # copy assembly helper scripts into conda env
plassembler download -d "$CONDA_PREFIX"/plassembler_db  # download plassembler_db
cd ../../bacteria_hybrid_assembly