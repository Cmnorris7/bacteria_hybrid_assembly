include { DNAAPLER_ALL }       from '../../modules/nf-core/dnaapler/all/main'
include { MEDAKA }             from '../../modules/nf-core/medaka/main'
include { GUNZIP_FASTA }       from '../../modules/local/gunzipfasta/main'
include { BWA_INDEX }          from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM }            from '../../modules/nf-core/bwa/mem/main'
include { POLYPOLISH_POLISH }  from '../../modules/nf-core/polypolish/polish/main'
include { PYPOLCA }            from '../../modules/local/pypolca/main'
include { SORT_FASTA }         from '../../modules/local/sortfasta/main'

workflow POLISHING {
    take:
    consensus
    reads

    main:
    consensus_ch = consensus.map { meta, assembly -> tuple(meta, assembly, []) }
    empty_db_ch = consensus.map { meta, assembly -> tuple([id: "${meta.id}_db"], []) }
    DNAAPLER_ALL(consensus_ch, empty_db_ch)

    long_reads_by_id_ch = reads.map { meta, short_reads, long_reads -> tuple(meta.id, long_reads) }
    medaka_input_ch = DNAAPLER_ALL.out.fasta
        .map { meta, assembly -> tuple(meta.id, meta, assembly) }
        .join(long_reads_by_id_ch)
        .map { id, meta, assembly, long_reads -> tuple(meta, long_reads, assembly) }
    MEDAKA(medaka_input_ch)
    GUNZIP_FASTA(MEDAKA.out.assembly)

    BWA_INDEX(GUNZIP_FASTA.out.fasta)

    short_reads_by_id_ch = reads.map { meta, short_reads, long_reads -> tuple(meta.id, short_reads) }
    bwa_jobs_ch = GUNZIP_FASTA.out.fasta
        .map { meta, assembly -> tuple(meta.id, meta, assembly) }
        .join(BWA_INDEX.out.index.map { meta, index -> tuple(meta.id, index) })
        .join(short_reads_by_id_ch)
        .flatMap { id, meta, assembly, index, short_reads ->
            short_reads.withIndex().collect { read, index_number ->
                def mate = index_number + 1
                def alignment_meta = meta + [sample_id: id, id: "${id}_R${mate}", single_end: true]
                tuple(alignment_meta, read, meta, index, assembly)
            }
        }

    bwa_inputs = bwa_jobs_ch.multiMap { alignment_meta, read, reference_meta, index, assembly ->
        reads: tuple(alignment_meta, read)
        index: tuple(reference_meta, index)
        fasta: tuple(reference_meta, assembly)
    }
    BWA_MEM(bwa_inputs.reads, bwa_inputs.index, bwa_inputs.fasta, false)

    sam_by_sample_ch = BWA_MEM.out.sam
        .map { meta, sam -> tuple(meta.sample_id, sam) }
        .groupTuple()
    polypolish_jobs_ch = GUNZIP_FASTA.out.fasta
        .map { meta, assembly -> tuple(meta.id, meta, assembly) }
        .join(sam_by_sample_ch)

    polypolish_inputs = polypolish_jobs_ch.multiMap { id, meta, assembly, sam_files ->
        assembly: tuple(meta, assembly)
        alignments: tuple(meta, sam_files)
    }
    POLYPOLISH_POLISH(polypolish_inputs.assembly, polypolish_inputs.alignments, false)

    paired_reads_by_id_ch = reads.map { meta, short_reads, long_reads -> tuple(meta.id, short_reads) }
    pypolca_input_ch = POLYPOLISH_POLISH.out.fasta
        .map { meta, assembly -> tuple(meta.id, meta, assembly) }
        .join(paired_reads_by_id_ch)
        .map { id, meta, assembly, short_reads -> tuple(meta, assembly, short_reads) }
    PYPOLCA(pypolca_input_ch)

    sort_script = file("${projectDir}/src/sort_fasta_add_length.py", checkIfExists: true)
    SORT_FASTA(PYPOLCA.out.assembly, sort_script)

    emit:
    assembly = SORT_FASTA.out.assembly
    dnaapler_summary = DNAAPLER_ALL.out.summary
    pypolca_log = PYPOLCA.out.log
}
