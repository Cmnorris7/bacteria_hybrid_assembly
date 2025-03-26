#!/bin/bash

#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name="bact_assembly_pipeline"
#SBATCH --partition=prod-compute
#SBATCH --output=slurm_pipeline_%j.log

# Bacterial Hybrid Assembly Pipeline - Master Script
# This script coordinates the execution of the complete assembly pipeline
module load slurm
module load nextflow

# Source the config file
export SCRIPT_DIR="$HOME/git/gitlab/bacteria_hybrid_assembly"
# add these to PATH ${HOME}/git/_github/Autocycler/scripts:${HOME}/git/_github/Autocycler/src:${HOME}/git/_github/Autocycler/target/release
export PATH="$PATH:${HOME}/git/_github/Autocycler/scripts:${HOME}/git/_github/Autocycler/src:${HOME}/git/_github/Autocycler/target/release"
echo "SCRIPT_DIR: $SCRIPT_DIR"
source "${SCRIPT_DIR}/config.sh"

# Exit on error
set -e

# Print commands before executing (for debugging)
set -x

# Set up log directory and file
log_file="$LOGS_DIR/run_bact_assembly.log"



# Log function
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

# Function to wait for a job to complete
wait_for_job() {
    local job_id=$1
    local step_name=$2
    
    log "Waiting for $step_name job (ID: $job_id) to complete..."
    
    # Wait for the job to finish
    while squeue -j "$job_id" -h -o "%t" | grep -q -E "PD|R|CF"; do
        sleep 60  # Check every minute
    done
    
    # Check if job completed successfully
    if sacct -j "$job_id" -o State -n | grep -q "COMPLETED"; then
        log "$step_name job completed successfully"
        return 0
    else
        log "ERROR: $step_name job failed. Check logs for details."
        return 1
    fi
}

# Create log directory if it doesn't exist
mkdir -p ${LOGS_DIR}

# Start pipeline
log "Starting bacterial hybrid assembly pipeline"

# double check that a file with 'nanopore' in the name exists, as well as 'R1' and 'R2' in current directory
if ! ls ${NANOPORE_PATTERN} &> /dev/null; then
    log "ERROR: No nanopore reads found in current directory. Must have 'nanopore' in the filename."
    exit 1
fi

if ! ls ${R1_PATTERN} &> /dev/null || ! ls ${R2_PATTERN} &> /dev/null; then
    log "ERROR: No Illumina reads found in current directory. Must have 'R1' and 'R2' in the filenames."
    exit 1
fi

# get sample name by splitting on _ and taking the first part
sample_name=$(ls ${NANOPORE_PATTERN} | cut -d'_' -f1)
file=$(ls ${NANOPORE_PATTERN})
# if sample_name is the same as the full file, just get basename of the file
if [ "$sample_name" == "$file" ]; then
    sample_name=$(basename "${file}" | cut -d'.' -f1)
fi

# if directory not the same as sample_name, make sample_name directory and move files to it
if [ "$(basename "$(pwd)")" != "$sample_name" ]; then
    log "Creating directory $sample_name and moving files into it"
    mkdir -p $sample_name
    mv *fastq* $sample_name
    cd $sample_name
    log_file="../$log_file"
fi

# Step 1: Read QC
log "Step 1: Running read quality control"
QC_JOB=$(sbatch --parsable --export=ALL ${SCRIPT_DIR}/readQC_sbatch.sh)
log "QC job submitted with ID: $QC_JOB"

# Wait for QC job to complete
if ! wait_for_job "$QC_JOB" "Read QC"; then
    log "Pipeline stopped due to failure in Read QC step"
    exit 1
fi

# Step 2: Autocycler assembly
log "Step 2: Running Autocycler assembly"
ASSEMBLY_JOB=$(sbatch --parsable --export=ALL ${SCRIPT_DIR}/autocycler_sbatch.sh)
log "Assembly job submitted with ID: $ASSEMBLY_JOB"

# Wait for assembly job to complete
if ! wait_for_job "$ASSEMBLY_JOB" "Autocycler assembly"; then
    log "Pipeline stopped due to failure in Autocycler assembly step"
    exit 1
fi

# Step 3: Polishing
log "Step 3: Running assembly polishing"
POLISH_JOB=$(sbatch --parsable --export=ALL ${SCRIPT_DIR}/polish_sbatch.sh)
log "Polishing job submitted with ID: $POLISH_JOB"

# Wait for polishing job to complete
if ! wait_for_job "$POLISH_JOB" "Assembly polishing"; then
    log "Pipeline stopped due to failure in Assembly polishing step"
    exit 1
fi

# Step 4: Assembly QC
log "Step 4: Running assembly quality assessment"
QC_ASSEMBLY_JOB=$(sbatch --parsable --export=ALL ${SCRIPT_DIR}/assemblyQC_sbatch.sh)
log "Assembly QC job submitted with ID: $QC_ASSEMBLY_JOB"

# Wait for assembly QC job to complete
if ! wait_for_job "$QC_ASSEMBLY_JOB" "Assembly quality assessment"; then
    log "Pipeline stopped due to failure in Assembly quality assessment step"
    exit 1
fi

# Pipeline complete
log "Pipeline completed successfully! Final assembly is in polish/*_final_assembly.fasta"
log "Assembly statistics are in nanopore_map_assembly_stats.tsv"

# Clean up intermediate files
rm ./fastp* 
mv autocycler.stderr ./autocycler_out
rm -r ./polish/tmp
mv ./slurm* ../${LOGS_DIR}
mv ../*log ../${LOGS_DIR}


