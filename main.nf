nextflow.enable.dsl = 2

include { INPUT_QC } from './subworkflows/local/input_qc'
include { ASSEMBLY_CONSENSUS } from './subworkflows/local/assembly_consensus'
include { POLISHING } from './subworkflows/local/polishing'
include { ASSEMBLY_QC } from './subworkflows/local/assembly_qc'

workflow {
    if (!params.input) {
        error "Missing required parameter: --input"
    }

    INPUT_QC(params.input)
    ASSEMBLY_CONSENSUS(INPUT_QC.out.reads)
    POLISHING(ASSEMBLY_CONSENSUS.out.assembly, INPUT_QC.out.reads)
    ASSEMBLY_QC(POLISHING.out.assembly, INPUT_QC.out.reads)
}
