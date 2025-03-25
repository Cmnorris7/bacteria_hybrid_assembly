#!/usr/bin/env bash

# This script is a wrapper for running a fully-automated Autocycler assembly.

# Usage:
#   autocycler_full.sh <read_fastq> <threads> <jobs>

# Copyright 2025 Ryan Wick (rrwick@gmail.com)
# Licensed under the GNU General Public License v3.
# See https://www.gnu.org/licenses/gpl-3.0.html.

# Ensure script exits on error.
set -e

# Get arguments.
reads=$1    # input reads FASTQ
threads=$2  # threads per job
jobs=$3     # number of simultaneous jobs

# Validate input parameters.
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    >&2 echo "Usage: $0 <read_fastq> <threads> <jobs>"
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    >&2 echo "Error: Input file '$reads' does not exist."
    exit 1
fi

# Function to wait for a job to complete and check its status
wait_for_job() {
    local job_id=$1
    local step_name=$2
    
    echo "Waiting for $step_name job (ID: $job_id) to complete..."
    
    # Wait for the job to finish
    while squeue -j "$job_id" -h -o "%t" | grep -q -E "PD|R|CF"; do
        sleep 60  # Check every minute
    done
    
    # Check if job completed successfully
    if sacct -j "$job_id" -o State -n | grep -q "COMPLETED"; then
        echo "$step_name job completed successfully"
        return 0
    else
        echo "ERROR: $step_name job failed. Check logs for details."
        return 1
    fi
}

# Function to wait for multiple jobs
wait_for_jobs() {
    local job_ids=("$@")
    local failed_jobs=0
    
    for job_id in "${job_ids[@]}"; do
        local state=$(sacct -j "$job_id" -n -o State | head -n 1 | awk '{print $1}')
        if [[ "$state" != "COMPLETED" ]]; then
            echo "Job $job_id failed with state: $state"
            failed_jobs=$((failed_jobs+1))
        fi
    done
    
    if [ $failed_jobs -gt 0 ]; then
        echo "Warning: $failed_jobs jobs failed."
        return 1
    else
        echo "All jobs completed successfully."
        return 0
    fi
}

genome_size=$(genome_size_raven.sh "$reads" "$threads")

# Step 1: subsample the long-read set into multiple files
echo "Step 1: Subsampling long reads"
autocycler subsample --reads "$reads" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file using SLURM
echo "Step 2: Assembling subsampled reads with multiple assemblers"
mkdir -p assemblies
job_ids=()

# Submit all assembly jobs
for i in 01 02 03 04; do
    echo "Submitting assembly jobs for sample $i"
    
    # Submit each assembler job and capture the job ID
    canu_id=$(sbatch --parsable --job-name=canu_"$i" --time=12:00:00 --mem=128000 --ntasks=1 --cpus-per-task=$threads \
              --wrap "canu.sh subsampled_reads/sample_$i.fastq assemblies/canu_$i $threads $genome_size")
    
    flye_id=$(sbatch --parsable --job-name=flye_"$i" --time=2:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
              --wrap "flye.sh subsampled_reads/sample_$i.fastq assemblies/flye_$i $threads $genome_size")
    
    miniasm_id=$(sbatch --parsable --job-name=miniasm_"$i" --time=1:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
                --wrap "miniasm.sh subsampled_reads/sample_$i.fastq assemblies/miniasm_$i $threads $genome_size")
    
    necat_id=$(sbatch --parsable --job-name=necat_"$i" --time=2:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
              --wrap "necat.sh subsampled_reads/sample_$i.fastq assemblies/necat_$i $threads $genome_size")
    
    nextdenovo_id=$(sbatch --parsable --job-name=nextdenovo_"$i" --time=2:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
                   --wrap "nextdenovo.sh subsampled_reads/sample_$i.fastq assemblies/nextdenovo_$i $threads $genome_size")
    
    raven_id=$(sbatch --parsable --job-name=raven_"$i" --time=1:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
              --wrap "raven.sh subsampled_reads/sample_$i.fastq assemblies/raven_$i $threads $genome_size")
    
    metamdbg_id=$(sbatch --parsable --job-name=metamdbg_"$i" --time=2:00:00 --mem=64000 --ntasks=1 --cpus-per-task=$threads \
                 --wrap "metamdbg.sh subsampled_reads/sample_$i.fastq assemblies/metamdbg_$i $threads $genome_size")
    
    # Add all job IDs to our array
    job_ids+=($canu_id $flye_id $miniasm_id $necat_id $nextdenovo_id $raven_id $metamdbg_id)
    
    echo "Submitted jobs for sample $i: $canu_id $flye_id $miniasm_id $necat_id $nextdenovo_id $raven_id $metamdbg_id"
done

# Create job dependency string with colon separators
job_dependency=$(IFS=:; echo "${job_ids[*]}")

# Submit a dependency job that will only run after all assembly jobs complete
echo "Submitting cleanup job dependent on all assembly jobs"
cleanup_job=$(sbatch --parsable --dependency=afterany:$job_dependency \
              --job-name="assembly_cleanup" --time=0:10:00 --mem=4G --ntasks=1 --wrap \
              "find assemblies/ -maxdepth 1 -type f -name \"*.fasta\" -empty -delete")

# Wait for the cleanup job to complete
wait_for_job "$cleanup_job" "Cleanup"

# Check status of all assembly jobs
echo "Checking status of all assembly jobs"
wait_for_jobs "${job_ids[@]}"
assembly_status=$?

if [ $assembly_status -ne 0 ]; then
    echo "Warning: Some assembly jobs failed. Pipeline will continue with available assemblies."
else
    echo "All assembly jobs completed successfully."
fi

# Optional step: remove the subsampled reads to save space
echo "Removing subsampled reads to save space"
rm -f subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
echo "Step 3: Compressing input assemblies into a unitig graph"
autocycler compress --max_contigs 40 -i assemblies -a autocycler_out 2>> autocycler.stderr

# Step 4: cluster the input contigs into putative genomic sequences
echo "Step 4: Clustering input contigs into putative genomic sequences"
autocycler cluster --max_contigs 40 -a autocycler_out 2>> autocycler.stderr

# Steps 5 and 6: trim and resolve each QC-pass cluster
echo "Steps 5 & 6: Trimming and resolving each QC-pass cluster"
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    cluster_name=$(basename "$c")
    echo "Processing $cluster_name"
    autocycler trim -c "$c" 2>> autocycler.stderr
    autocycler resolve -c "$c" 2>> autocycler.stderr
done

# Step 7: combine resolved clusters into a final assembly
echo "Step 7: Combining resolved clusters into final assembly"
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr

echo "Assembly pipeline completed successfully"