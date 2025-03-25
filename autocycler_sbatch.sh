#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --job-name="bact_assembly_autocycler"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --export=NONE
#SBATCH --output=slurm_autocycler_%j.log


# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "FILTERED_READS_DIR: $FILTERED_READS_DIR"

. ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${AUTOCYCLER_ENV}
module load parallel

# Get the nanopore filtered reads
nanopore=${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz

${SCRIPT_DIR}/autocycler_assembly_pipeline.sh ${nanopore} 16 4