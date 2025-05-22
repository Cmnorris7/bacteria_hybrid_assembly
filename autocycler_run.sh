#!/bin/bash

# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "FILTERED_READS_DIR: $FILTERED_READS_DIR"

. ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${AUTOCYCLER_ENV}

# Get the nanopore filtered reads
nanopore=${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz

# Calculate threads divisible by 4 (rounded down)
THREADS_DIV4=$(( $THREADS / 4 ))
ASSEMBLER_THREADS=4  # Fixed number of assemblers to use (job count)

# Log thread allocation
echo "Total threads: $THREADS"
echo "Threads per job: $THREADS_DIV4"
echo "Number of parallel jobs: $ASSEMBLER_THREADS"

${AUTOCYCLER_DIR}/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh ${nanopore} $THREADS_DIV4 $ASSEMBLER_THREADS