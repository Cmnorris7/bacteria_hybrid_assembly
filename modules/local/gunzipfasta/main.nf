process GUNZIP_FASTA {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/pypolca@sha256:82af8f2cba0f8c9aa8c1d221e0ac16bba423ba010ebba5e98840d08405f63519'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.medaka"
    """
    gzip -dc ${assembly} > ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.medaka"
    """
    touch ${prefix}.fasta
    """
}
