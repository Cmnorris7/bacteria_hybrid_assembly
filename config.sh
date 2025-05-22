#!/bin/bash

# Script directory - the directory where all pipeline scripts are located
export SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
export AUTOCYCLER_DIR="$(readlink -f "${SCRIPT_DIR}/../Autocycler")"

# Pipeline directories
export POLISH_DIR="polish"
export FILTERED_READS_DIR="filtered_reads"
export ASSEMBLY_QC_DIR="assembly_qc"
export AUTOCYCLER_OUT_DIR="autocycler_out"
export TMP_DIR="tmp"
export LOGS_DIR="logs"

# Conda environments
export CONDA_BASE=$(conda info --base)
export AUTOCYCLER_ENV="${CONDA_BASE}/envs/autocycler"
export QC_ENV="${CONDA_BASE}/envs/qc_polish_hybrid"

# Custom scripts
export MAP_STATS_SCRIPT="${SCRIPT_DIR}/assemblyQC_coverage_stats.sh"
export SORT_FASTA_SCRIPT="${SCRIPT_DIR}/sort_fasta_add_length.py"

# Workflow parameters
export MEDAKA_MODEL="r1041_e82_400bps_hac_v4.3.0"
export THREADS=$(nproc --all) # will default to all available threads
# export THREADS=24
export SHORT_READ_LENGTH=90

# File patterns
export NANOPORE_PATTERN="*nanopore*fastq.gz"
export R1_PATTERN="*R1*fastq.gz"
export R2_PATTERN="*R2*fastq.gz"
export FINAL_ASSEMBLY_PATTERN="*final*"

