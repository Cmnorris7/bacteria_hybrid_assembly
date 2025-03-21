#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --job-name="assemblyQC"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --export=NONE
#SBATCH --output=slurm_assemblyQC_%j.log

# Source the config file
export SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
source "${SCRIPT_DIR}/config.sh"

module load bwa/0.7.17
module load minimap2-2.24-gcc-9.2.0-xnxd4
module load samtools-1.12-gcc-9.2.0-m7qh35t 

# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "POLISH_DIR: $POLISH_DIR"
echo "ASSEMBLY_QC_DIR: $ASSEMBLY_QC_DIR"

# Get the final assembly file
final=${POLISH_DIR}/${SAMPLE_NAME}_final_assembly.fasta

# Get the nanopore file from filtered_reads directory
nanopore=${FILTERED_READS_DIR}/${SAMPLE_NAME}_nanopore_filtered.fastq.gz

# Create output directory
mkdir -p ${ASSEMBLY_QC_DIR}

# index final with both bwa and minimap2
samtools faidx ${final}

# long read mapping to assembly
minimap2 -t ${THREADS} -ax map-ont ${final} ${nanopore} | samtools view -bS - | samtools sort -o ${ASSEMBLY_QC_DIR}/nanopore.sorted.bam
${MAP_STATS_SCRIPT} ${ASSEMBLY_QC_DIR}/nanopore.sorted.bam ${ASSEMBLY_QC_DIR}/nanopore_map_assembly_stats.tsv

# Copy final stats to main directory for easy access
cp ${ASSEMBLY_QC_DIR}/nanopore_map_assembly_stats.tsv ./nanopore_map_assembly_stats.tsv

# if curious, can also align & assess short reads
# uncomment the following lines

# R1=${FILTERED_READS_DIR}/${SAMPLE_NAME}_R1_filtered.fastq.gz
# R2=${FILTERED_READS_DIR}/${SAMPLE_NAME}_R2_filtered.fastq.gz
# bwa index ${final}
# bwa mem -t ${THREADS} ${final} ${R1} ${R2} | samtools view -bS - | samtools sort -o ${ASSEMBLY_QC_DIR}/illumina.sorted.bam
# ${MAP_STATS_SCRIPT} ${ASSEMBLY_QC_DIR}/illumina.sorted.bam ${ASSEMBLY_QC_DIR}/illumina_map_assembly_stats.tsv

