#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --job-name="bact_assembly_polish"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --output=slurm_polish_%j.log
#SBATCH --export=NONE

# Workflow summary
# 1. Copy the consensus assembly from autocycler_out to the polish directory
# 2. Rotate the assembly
# 3. Run medaka consensus on the long reads
# 4. Run polypolish using the short reads
# 5. Run pypolca using the short reads
# 6. Fix the fasta naming


######### Setup #########
. ${CONDA_BASE}/etc/profile.d/conda.sh

# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "FILTERED_READS_DIR: $FILTERED_READS_DIR"
echo "POLISH_DIR: $POLISH_DIR"
echo "AUTOCYCLER_OUT_DIR: $AUTOCYCLER_OUT_DIR"

dir_name=$(basename "$(pwd)")
mkdir -p ${POLISH_DIR}
cd ${POLISH_DIR}
mkdir -p ${TMP_DIR}

######### Copy the consensus assembly from autocycler_out #########
cp ../${AUTOCYCLER_OUT_DIR}/consensus_assembly.fasta ./${TMP_DIR}/draft.fasta

######### Rotate the assembly #########
singularity exec ${DNAAPLER_CONTAINER} dnaapler all --autocomplete mystery --seed_value 13 -i ./${TMP_DIR}/draft.fasta -o ./dnaapler/pre-polish/ -t ${THREADS}
rm ./${TMP_DIR}/draft.fasta
cp ./dnaapler/pre-polish/dnaapler_reoriented.fasta ./${TMP_DIR}/draft_dnaapler.fasta

######### Run medaka long read polishing#########
# may need to change the model depending on the data
singularity exec ${MEDAKA_CONTAINER} medaka_consensus -i ../${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz -d ./${TMP_DIR}/draft_dnaapler.fasta -o medaka -t ${THREADS} -m ${MEDAKA_MODEL} --bacteria
cp ./medaka/consensus.fasta ./${TMP_DIR}/draft_medaka.fasta

######### Run polypolish short read polishing#########
conda activate ${POLYPOLISH_ENV}

bwa index ./${TMP_DIR}/draft_medaka.fasta
bwa mem -t ${THREADS} -a ./${TMP_DIR}/draft_medaka.fasta ../${FILTERED_READS_DIR}/R1_filtered.fastq.gz > alignments_1.sam
bwa mem -t ${THREADS} -a ./${TMP_DIR}/draft_medaka.fasta ../${FILTERED_READS_DIR}/R2_filtered.fastq.gz > alignments_2.sam
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish ./${TMP_DIR}/draft_medaka.fasta filtered_1.sam filtered_2.sam > ./${TMP_DIR}/draft_polypolish.fasta
rm *.amb *.ann *.bwt *.pac *.sa *.sam

######### Run pypolca short read polishing#########
conda activate ${PYPOLCA_ENV}

pypolca run -a ./${TMP_DIR}/draft_polypolish.fasta -1 ../${FILTERED_READS_DIR}/R1_filtered.fastq.gz -2 ../${FILTERED_READS_DIR}/R2_filtered.fastq.gz -t ${THREADS} -o pypolca --careful
cp ./pypolca/corrected.fasta ./${TMP_DIR}/draft_pypolca.fasta

######### Fix fasta naming#########
# run the script to fix the fasta naming and add lengths
python3 ${SORT_FASTA_SCRIPT} ./${TMP_DIR}/draft_pypolca.fasta ./final_assembly.fasta