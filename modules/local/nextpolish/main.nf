process NEXTPOLISH {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/nextpolish@sha256:d2c4525eb4aa4610fd42a18d3eafc1005e57f4469bfbbd6458e6bb9c7e1cc968'

    input:
    tuple val(meta), path(assembly), path(reads)
    val read_type

    output:
    tuple val(meta), path('*.fasta.gz'), emit: assembly
    tuple val(meta), path('*.nextpolish.log'), emit: log
    tuple val("${task.process}"), val('nextpolish'), eval('nextPolish -v 2>&1 | awk \'{print $2}\''), emit: versions_nextpolish, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def settings = [
        ont:  [section: 'lgs',  preset: 'map-ont'],
        clr:  [section: 'lgs',  preset: 'map-pb'],
        hifi: [section: 'hifi', preset: 'map-hifi']
    ][read_type]
    if (!settings) {
        error "Unsupported NextPolish read type '${read_type}'; expected ont, clr, or hifi"
    }
    """
    set -o pipefail
    printf '%s\n' "\$(readlink -f ${reads})" > input.fofn

    printf '%s\n' \
        '[General]' \
        'job_type = local' \
        'job_prefix = nextPolish' \
        'task = best' \
        'rewrite = yes' \
        'deltmp = yes' \
        'rerun = 3' \
        'parallel_jobs = 1' \
        'multithread_jobs = ${task.cpus}' \
        'genome = ${assembly}' \
        'genome_size = auto' \
        'workdir = nextpolish' \
        'polish_options = -p ${task.cpus}' \
        '' \
        '[${settings.section}_option]' \
        '${settings.section}_fofn = input.fofn' \
        '${settings.section}_options = -min_read_len 1k -max_depth 100' \
        '${settings.section}_minimap2_options = -x ${settings.preset} -t ${task.cpus}' > nextpolish_run.cfg

    nextPolish nextpolish_run.cfg ${args} 2>&1 | tee ${prefix}.nextpolish.log
    test -s nextpolish/genome.nextpolish.fasta
    gzip -c nextpolish/genome.nextpolish.fasta > ${prefix}.fasta.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip -c > ${prefix}.fasta.gz
    touch ${prefix}.nextpolish.log
    """
}