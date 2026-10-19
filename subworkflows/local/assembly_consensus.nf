include { AUTOCYCLER_SUBSAMPLE } from '../../modules/nf-core/autocycler/subsample/main'
include { FLYE }                 from '../../modules/nf-core/flye/main'
include { MINIMAP2_ALIGN }       from '../../modules/nf-core/minimap2/align/main'
include { MINIASM }              from '../../modules/nf-core/miniasm/main'
include { MINIPOLISH }           from '../../modules/local/minipolish/main'
include { RAVEN }                from '../../modules/nf-core/raven/main'
include { METAMDBG_ASM }         from '../../modules/nf-core/metamdbg/asm/main'
include { CANU }                 from '../../modules/nf-core/canu/main'
include { NECAT }                from '../../modules/local/necat/main'
include { NEXTDENOVO }           from '../../modules/local/nextdenovo/main'
include { NEXTPOLISH }           from '../../modules/local/nextpolish/main'
include { AUTOCYCLER_COMPRESS }  from '../../modules/nf-core/autocycler/compress/main'
include { AUTOCYCLER_CLUSTER }   from '../../modules/nf-core/autocycler/cluster/main'
include { AUTOCYCLER_TRIM }      from '../../modules/nf-core/autocycler/trim/main'
include { AUTOCYCLER_RESOLVE }   from '../../modules/nf-core/autocycler/resolve/main'
include { AUTOCYCLER_COMBINE }   from '../../modules/nf-core/autocycler/combine/main'

workflow ASSEMBLY_CONSENSUS {
    take:
    reads

    main:
    enabled_assemblers = params.assemblers.tokenize(',')*.trim() as Set
    supported_assemblers = ['flye', 'miniasm', 'raven', 'metamdbg', 'canu', 'necat', 'nextdenovo'] as Set
    unknown_assemblers = enabled_assemblers - supported_assemblers
    if (unknown_assemblers) {
        error "Unsupported assembler(s): ${unknown_assemblers.sort().join(', ')}"
    }
    if (enabled_assemblers.size() < 3) {
        error 'Autocycler requires at least three enabled assemblers'
    }

    long_reads_ch = reads.map { meta, short_reads, long_reads -> tuple(meta, long_reads) }
    AUTOCYCLER_SUBSAMPLE(long_reads_ch, params.genome_size)

    subsamples_ch = AUTOCYCLER_SUBSAMPLE.out.subsampled_reads.flatMap { meta, subset_reads ->
        def subsets = subset_reads instanceof List ? subset_reads : [subset_reads]
        subsets.withIndex().collect { subset, index ->
            def subset_number = index + 1
            def subset_meta = meta + [
                sample_id: meta.id,
                subset: subset_number,
                id: "${meta.id}_subsample_${String.format('%02d', subset_number)}"
            ]
            tuple(subset_meta, subset)
        }
    }

    flye_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('flye') }
    FLYE(flye_reads_ch, params.long_read_mode == 'nano-hq' ? '--nano-hq' : '--nano-raw')

    miniasm_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('miniasm') }
    miniasm_reference_ch = miniasm_reads_ch.map { meta, subset -> tuple(meta, []) }
    MINIMAP2_ALIGN(miniasm_reads_ch, miniasm_reference_ch, false, '', false, false)
    miniasm_input_ch = miniasm_reads_ch
        .join(MINIMAP2_ALIGN.out.paf)
        .map { meta, subset, paf -> tuple(meta, subset, paf) }
    MINIASM(miniasm_input_ch)
    minipolish_input_ch = miniasm_reads_ch
        .join(MINIASM.out.gfa)
        .map { meta, subset, gfa -> tuple(meta, subset, gfa) }
    MINIPOLISH(minipolish_input_ch, params.long_read_mode == 'nano-hq' ? 'lr:hq' : 'map-ont')

    raven_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('raven') }
    RAVEN(raven_reads_ch)

    metamdbg_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('metamdbg') }
    METAMDBG_ASM(metamdbg_reads_ch, 'ont')

    canu_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('canu') }
    CANU(canu_reads_ch, '-nanopore', params.genome_size)

    necat_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('necat') }
    NECAT(necat_reads_ch, params.genome_size)

    nextdenovo_reads_ch = subsamples_ch.filter { meta, subset -> enabled_assemblers.contains('nextdenovo') }
    NEXTDENOVO(nextdenovo_reads_ch, params.genome_size, 'ont')
    nextpolish_input_ch = NEXTDENOVO.out.assembly
        .join(nextdenovo_reads_ch)
        .map { meta, assembly, subset -> tuple(meta, assembly, subset) }
    NEXTPOLISH(nextpolish_input_ch, 'ont')

    candidate_assemblies_ch = FLYE.out.fasta
        .mix(MINIPOLISH.out.assembly)
        .mix(RAVEN.out.fasta)
        .mix(METAMDBG_ASM.out.contigs)
        .mix(CANU.out.assembly)
        .mix(NECAT.out.assembly)
        .mix(NEXTPOLISH.out.assembly)
        .map { meta, assembly -> tuple([id: meta.sample_id], assembly) }
        .groupTuple()

    AUTOCYCLER_COMPRESS(candidate_assemblies_ch)
    AUTOCYCLER_CLUSTER(AUTOCYCLER_COMPRESS.out.gfa)

    cluster_graphs_ch = AUTOCYCLER_CLUSTER.out.clusters.flatMap { meta, graphs ->
        def graph_list = graphs instanceof List ? graphs : [graphs]
        graph_list.collect { graph ->
            tuple(meta + [sample_id: meta.id, id: "${meta.id}_${graph.simpleName}"], graph)
        }
    }
    AUTOCYCLER_TRIM(cluster_graphs_ch)
    AUTOCYCLER_RESOLVE(AUTOCYCLER_TRIM.out.gfa)

    resolved_clusters_ch = AUTOCYCLER_RESOLVE.out.resolved
        .map { meta, graph -> tuple([id: meta.sample_id], graph) }
        .groupTuple()
    AUTOCYCLER_COMBINE(resolved_clusters_ch)

    emit:
    assembly = AUTOCYCLER_COMBINE.out.fasta
    graph = AUTOCYCLER_COMBINE.out.gfa
    stats = AUTOCYCLER_COMBINE.out.stats
}
