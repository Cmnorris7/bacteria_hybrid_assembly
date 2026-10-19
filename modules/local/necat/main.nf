process NECAT {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/necat@sha256:6a13431f8c8768ccf30830206af534a413baa6f5a5c6a917a7729c70729d53cd'

    input:
    tuple val(meta), path(reads)
    val genome_size

    output:
    tuple val(meta), path('*.fasta.gz'), emit: assembly
    tuple val(meta), path('*.necat.log'), emit: log
    tuple val("${task.process}"), val('necat'), val('0.0.1_update20200803'), emit: versions_necat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\n' "\$(readlink -f ${reads})" > read_list.txt

    printf '%s\n' \
        'PROJECT=necat' \
        'ONT_READ_LIST=read_list.txt' \
        'GENOME_SIZE=${genome_size}' \
        'THREADS=${task.cpus}' \
        'MIN_READ_LENGTH=3000' \
        'PREP_OUTPUT_COVERAGE=40' \
        'OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000' \
        'OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000' \
        'CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0' \
        'CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0' \
        'TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400' \
        'ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400' \
        'NUM_ITER=2' \
        'CNS_OUTPUT_COVERAGE=30' \
        'CLEANUP=1' \
        'USE_GRID=false' \
        'GRID_NODE=0' \
        'GRID_OPTIONS=' \
        'SMALL_MEMORY=0' \
        'FSA_OL_FILTER_OPTIONS=' \
        'FSA_ASSEMBLE_OPTIONS=' \
        'FSA_CTG_BRIDGE_OPTIONS=' \
        'POLISH_CONTIGS=true' > config.txt

    necat bridge config.txt ${args} 2>&1 | tee ${prefix}.necat.log
    gzip -c necat/6-bridge_contigs/polished_contigs.fasta > ${prefix}.fasta.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip -c > ${prefix}.fasta.gz
    touch ${prefix}.necat.log
    """
}