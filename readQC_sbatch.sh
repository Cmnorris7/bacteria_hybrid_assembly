#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --job-name="bact_assembly_readQC"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --export=NONE
#SBATCH --output=slurm_readQC_%j.log


#################### Env prep ####################
. ${CONDA_BASE}/etc/profile.d/conda.sh
mkdir -p ${FILTERED_READS_DIR}

# Find input files directly using patterns
NANOPORE_FILE=$(ls ${NANOPORE_PATTERN} 2>/dev/null | head -n1)
R1_FILE=$(ls ${R1_PATTERN} 2>/dev/null | head -n1)
R2_FILE=$(ls ${R2_PATTERN} 2>/dev/null | head -n1)

# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "FILTERED_READS_DIR: $FILTERED_READS_DIR"
echo "THREADS: $THREADS"
echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "Input files:"
echo "  Nanopore: $NANOPORE_FILE"
echo "  R1: $R1_FILE"
echo "  R2: $R2_FILE"

#################### Fastp for Illumina pe reads ####################
conda activate ${CONDA_BASE}/envs/fastp
echo "Using fastp to filter Illumina reads.."
fastp -i ${R1_FILE} -I ${R2_FILE} \
      -o ./${FILTERED_READS_DIR}/R1_filtered.fastq.gz \
      -O ./${FILTERED_READS_DIR}/R2_filtered.fastq.gz \
      --length_required 140 --thread $((THREADS-2))
echo "fastp complete."

#################### Filtlong for nanopore reads ####################
echo "Using filtlong to filter Nanopore reads.."
# if illumina reads are good, use this command
# ~/git/_github/Filtlong/bin/filtlong -1 ./${FILTERED_READS_DIR}/R1_filtered.fastq.gz \
#   -2 ./${FILTERED_READS_DIR}/R2_filtered.fastq.gz \
#   --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 500 \
#   ${NANOPORE_FILE} | gzip > ./${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz

# if illumina reads are bad, use this command
## defaulting to this command to be conservative
~/git/_github/Filtlong/bin/filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 \
  ${NANOPORE_FILE} | gzip > ./${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz
echo "filtlong complete."

#################### FASTQC plots ####################
echo "Making fastqc plots for filtered reads.."
conda activate ${CONDA_BASE}/envs/fastqc
mkdir -p ${FILTERED_READS_DIR}/fastqc
fastqc ./${FILTERED_READS_DIR}/*fastq.gz -o ./${FILTERED_READS_DIR}/fastqc --threads $((THREADS-2))

#################### Clean up ####################
mkdir -p raw_reads
mv ./*fastq.gz ./raw_reads
echo "Filtered reads moved to ${FILTERED_READS_DIR} directory. Raw reads moved to raw_reads directory."

echo "Hybrid QC complete."