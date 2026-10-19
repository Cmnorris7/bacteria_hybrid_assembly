process SORT_FASTA {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/pypolca@sha256:82af8f2cba0f8c9aa8c1d221e0ac16bba423ba010ebba5e98840d08405f63519'

    input:
    tuple val(meta), path(assembly)
    path sort_script

    output:
    tuple val(meta), path('*_final_assembly.fasta'), emit: assembly

    script:
    """
    python ${sort_script} ${assembly} ${meta.id}_final_assembly.fasta
    """

    stub:
    """
    touch ${meta.id}_final_assembly.fasta
    """
}
