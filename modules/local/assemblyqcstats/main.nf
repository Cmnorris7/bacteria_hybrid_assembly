process ASSEMBLY_QC_STATS {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    input:
    tuple val(meta), path(bam), path(bai)
    path stats_script

    output:
    tuple val(meta), path('*.mapping_stats.tsv'), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bash ${stats_script} ${bam} ${prefix}.mapping_stats.tsv ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Metric\tValue\nTotal_Reads\t0\nMapped_Reads\t0\nMapping_Rate_Percent\t0.00\n' > ${prefix}.mapping_stats.tsv
    """
}