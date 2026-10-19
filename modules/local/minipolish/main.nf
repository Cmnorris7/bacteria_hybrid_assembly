process MINIPOLISH {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/minipolish@sha256:d1f78f37d189627b1d7d9b844c3f3765be422a224bc946f03f27dd97912c1c32'

    input:
    tuple val(meta), path(reads), path(gfa)
    val minimap2_preset

    output:
    tuple val(meta), path('*.fasta.gz'), emit: assembly
    tuple val(meta), path('*.gfa.gz'), emit: gfa
    tuple val("${task.process}"), val('minipolish'), eval("minipolish --version 2>&1 | sed -n 's/^Minipolish v//p'"), emit: versions_minipolish, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def minDepthRel = task.ext.min_depth_rel ?: 0.0
    def prefix = task.ext.prefix ?: "${meta.id}.minipolish"
    """
    set -o pipefail
    gzip -dc ${gfa} > input.gfa

    minipolish \
        --threads ${task.cpus} \
        --minimap2-preset ${minimap2_preset} \
        ${args} \
        ${reads} \
        input.gfa > ${prefix}.gfa

    awk -F '\t' -v min_depth_rel=${minDepthRel} '
        \$1 == "S" {
            count++
            names[count] = \$2
            sequences[count] = \$3
            depths[count] = ""
            for (i = 4; i <= NF; i++) {
                if (\$i ~ /^dp:f:/ || \$i ~ /^rd:i:/) {
                    split(\$i, depth, ":")
                    depths[count] = depth[3]
                    break
                }
            }
            if (length(\$3) > longest_length) {
                longest_length = length(\$3)
                longest_depth = depths[count]
            }
        }
        END {
            filter = min_depth_rel > 0 && longest_depth != ""
            if (filter) {
                for (i = 1; i <= count; i++) {
                    if (depths[i] == "") filter = 0
                }
            }
            for (i = 1; i <= count; i++) {
                if (filter && depths[i] < longest_depth * min_depth_rel) continue
                header = names[i]
                if (names[i] ~ /c\$/) header = header " circular=true"
                if (depths[i] != "") header = header " depth=" depths[i]
                print ">" header
                print sequences[i]
            }
        }
    ' ${prefix}.gfa > ${prefix}.fasta

    test -s ${prefix}.fasta
    gzip -n ${prefix}.gfa
    gzip -n ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.minipolish"
    """
    echo '' | gzip -c > ${prefix}.fasta.gz
    echo '' | gzip -c > ${prefix}.gfa.gz
    """
}