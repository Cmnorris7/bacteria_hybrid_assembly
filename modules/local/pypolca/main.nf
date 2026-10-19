process PYPOLCA {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/pypolca@sha256:82af8f2cba0f8c9aa8c1d221e0ac16bba423ba010ebba5e98840d08405f63519'

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path('*.pypolca.fasta'), emit: assembly
    tuple val(meta), path('*.pypolca.log'), emit: log
    tuple val("${task.process}"), val('pypolca'), eval('pypolca --version 2>&1'), emit: versions_pypolca, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pypolca run \
        -a ${assembly} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -t ${task.cpus} \
        -o pypolca \
        --careful \
        ${args} 2>&1 | tee ${prefix}.pypolca.log

    cp pypolca/pypolca_corrected.fasta ${prefix}.pypolca.fasta
    """

            stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pypolca.fasta
    touch ${prefix}.pypolca.log
    """
}
