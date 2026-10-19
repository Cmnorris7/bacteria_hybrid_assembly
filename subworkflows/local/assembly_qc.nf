include { MINIMAP2_ALIGN as MINIMAP2_QC } from '../../modules/nf-core/minimap2/align/main'
include { ASSEMBLY_QC_STATS }                  from '../../modules/local/assemblyqcstats/main'

workflow ASSEMBLY_QC {
    take:
    assembly
    reads

    main:
    long_reads_by_id_ch = reads.map { meta, short_reads, long_reads -> tuple(meta.id, long_reads) }
    qc_jobs_ch = assembly
        .map { meta, fasta -> tuple(meta.id, meta, fasta) }
        .join(long_reads_by_id_ch)

    mapping_inputs = qc_jobs_ch.multiMap { id, meta, fasta, long_reads ->
        reads: tuple(meta, long_reads)
        reference: tuple(meta, fasta)
    }
    MINIMAP2_QC(mapping_inputs.reads, mapping_inputs.reference, true, 'bai', false, false)

    stats_script = file("${projectDir}/src/assemblyQC_coverage_stats.sh", checkIfExists: true)
    stats_input_ch = MINIMAP2_QC.out.bam
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(MINIMAP2_QC.out.index.map { meta, bai -> tuple(meta.id, bai) })
        .map { id, meta, bam, bai -> tuple(meta, bam, bai) }
    ASSEMBLY_QC_STATS(stats_input_ch, stats_script)

    emit:
    bam = MINIMAP2_QC.out.bam
    bai = MINIMAP2_QC.out.index
    stats = ASSEMBLY_QC_STATS.out.stats
}