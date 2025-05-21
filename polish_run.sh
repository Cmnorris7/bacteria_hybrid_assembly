#!/bin/bash

# Workflow summary
# 1. Copy the consensus assembly from autocycler_out to the polish directory
# 2. Rotate the assembly
# 3. Run medaka consensus on the long reads
# 4. Run polypolish using the short reads
# 5. Run pypolca using the short reads
# 6. Fix the fasta naming


######### Setup #########
. ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${QC_ENV}

# Display important environment variables for debugging
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "FILTERED_READS_DIR: $FILTERED_READS_DIR"
echo "POLISH_DIR: $POLISH_DIR"
echo "AUTOCYCLER_OUT_DIR: $AUTOCYCLER_OUT_DIR"

dir_name=$(basename "$(pwd)")
mkdir -p ${POLISH_DIR}
# cd ${POLISH_DIR}
WORKDIR=${POLISH_DIR}/${TMP_DIR}
mkdir -p ${WORKDIR}

######### Copy the consensus assembly from autocycler_out #########
cp ./${AUTOCYCLER_OUT_DIR}/consensus_assembly.fasta ./${WORKDIR}/draft.fasta

######### Rotate the assembly #########
dnaapler all --autocomplete mystery --seed_value 13 -i ./${WORKDIR}/draft.fasta -o ${POLISH_DIR}/dnaapler/pre-polish/ -t $((THREADS - 2))}
rm ${WORKDIR}/draft.fasta
cp ${POLISH_DIR}/dnaapler/pre-polish/dnaapler_reoriented.fasta ./${WORKDIR}/draft_dnaapler.fasta

######### Run medaka long read polishing#########
# may need to change the model depending on the data
medaka_consensus -i ${FILTERED_READS_DIR}/nanopore_filtered.fastq.gz -d ${WORKDIR}/draft_dnaapler.fasta -o ${POLISH_DIR}/medaka -t $((THREADS - 2))} -m ${MEDAKA_MODEL} --bacteria
cp ${POLISH_DIR}/medaka/consensus.fasta ./${WORKDIR}/draft_medaka.fasta

######### Run polypolish short read polishing#########
bwa index ./${WORKDIR}/draft_medaka.fasta
bwa mem -t $((THREADS - 2))} -a ${WORKDIR}/draft_medaka.fasta ${FILTERED_READS_DIR}/R1_filtered.fastq.gz > ${WORKDIR}/alignments_1.sam
bwa mem -t $((THREADS - 2))} -a ${WORKDIR}/draft_medaka.fasta ${FILTERED_READS_DIR}/R2_filtered.fastq.gz > ${WORKDIR}/alignments_2.sam
polypolish filter --in1 ${WORKDIR}/alignments_1.sam --in2 ${WORKDIR}/alignments_2.sam --out1 ${WORKDIR}/filtered_1.sam --out2 ${WORKDIR}/filtered_2.sam
polypolish polish ${WORKDIR}/draft_medaka.fasta ${WORKDIR}/filtered_1.sam ${WORKDIR}/filtered_2.sam > ${WORKDIR}/draft_polypolish.fasta

######### Run pypolca short read polishing#########
pypolca run -a ./${WORKDIR}/draft_polypolish.fasta -1 ${FILTERED_READS_DIR}/R1_filtered.fastq.gz -2 ${FILTERED_READS_DIR}/R2_filtered.fastq.gz -t $((THREADS - 2))} -o ${POLISH_DIR}/pypolca --careful
cp ${POLISH_DIR}/pypolca/pypolca_corrected.fasta ./${WORKDIR}/draft_pypolca.fasta

######### Fix fasta naming#########
# run the script to fix the fasta naming and add lengths
python3 ${SORT_FASTA_SCRIPT} ${WORKDIR}/draft_pypolca.fasta ${POLISH_DIR}/final_assembly.fasta