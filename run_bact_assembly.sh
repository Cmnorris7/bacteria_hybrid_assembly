#!/bin/bash

# Bacterial Hybrid Assembly Pipeline - Master Script

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
log "Starting Read QC process locally"

# Execute the script directly
bash ${SCRIPT_DIR}/readQC_run.sh

# Check if the command completed successfully
if [ $? -eq 0 ]; then
    log "Read QC completed successfully"
else
    log "Pipeline stopped due to failure in Read QC step"
    exit 1
fi

# Step 2: Autocycler assembly
log "Step 2: Running Autocycler assembly"
log "Starting Autocycler assembly process locally"

# Execute the script directly
bash ${SCRIPT_DIR}/autocycler_run.sh

# Check if the command completed successfully
if [ $? -eq 0 ]; then
    log "Autocycler assembly completed successfully"
else
    log "Pipeline stopped due to failure in Autocycler assembly step"
    exit 1
fi

# Step 3: Polishing
log "Step 3: Running assembly polishing"
log "Starting assembly polishing process locally"

# Execute the script directly
bash ${SCRIPT_DIR}/polish_run.sh

# Check if the command completed successfully
if [ $? -eq 0 ]; then
    log "Assembly polishing completed successfully"
else
    log "Pipeline stopped due to failure in Assembly polishing step"
    exit 1
fi

# Step 4: Assembly QC
log "Step 4: Running assembly quality assessment"
log "Starting assembly quality assessment process locally"

# Execute the script directly
bash ${SCRIPT_DIR}/assemblyQC_run.sh

# Check if the command completed successfully
if [ $? -eq 0 ]; then
    log "Assembly quality assessment completed successfully"
else
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