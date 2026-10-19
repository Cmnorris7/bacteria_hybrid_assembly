process NEXTDENOVO {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/nextdenovo@sha256:32180f4c243488f70e1e11393d86fa99f4f5ea0b5a0b7d8c770b3f9a9dbe191a'

    input:
    tuple val(meta), path(reads)
    val genome_size
    val read_type

    output:
    tuple val(meta), path('*.nextdenovo.fasta'), emit: assembly
    tuple val(meta), path('*.nextdenovo.log'), emit: log
    tuple val("${task.process}"), val('nextdenovo'), eval('nextDenovo -v 2>&1 | awk \'{print $2}\''), emit: versions_nextdenovo, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail
    printf '%s\n' "\$(readlink -f ${reads})" > input.fofn

    printf '%s\n' \
        '[General]' \
        'job_type = local' \
        'job_prefix = nextDenovo' \
        'task = all' \
        'rewrite = yes' \
        'deltmp = yes' \
        'parallel_jobs = 1' \
        'input_type = raw' \
        'read_type = ${read_type}' \
        'input_fofn = input.fofn' \
        'workdir = nextdenovo' \
        '' \
        '[correct_option]' \
        'read_cutoff = 1k' \
        'genome_size = ${genome_size}' \
        'sort_options = -m 20g -t ${task.cpus}' \
        'minimap2_options_raw = -t ${task.cpus}' \
        'pa_correction = 1' \
        'correction_options = -p ${task.cpus}' \
        '' \
        '[assemble_option]' \
        'minimap2_options_cns = -t ${task.cpus}' \
        'nextgraph_options = -a 1' > nextdenovo_run.cfg

    nextDenovo nextdenovo_run.cfg ${args} 2>&1 | tee ${prefix}.nextdenovo.log
    test -s nextdenovo/03.ctg_graph/nd.asm.fasta
    cp nextdenovo/03.ctg_graph/nd.asm.fasta ${prefix}.nextdenovo.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nextdenovo.fasta
    touch ${prefix}.nextdenovo.log
    """
}