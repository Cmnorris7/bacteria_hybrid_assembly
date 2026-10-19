include { FASTP }    from '../../modules/nf-core/fastp/main'
include { FILTLONG } from '../../modules/nf-core/filtlong/main'
include { FASTQC }   from '../../modules/nf-core/fastqc/main'

workflow INPUT_QC {
    take:
    samplesheet

    main:
    samples_ch = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.long_reads || !row.short_reads_1 || !row.short_reads_2) {
                error "Samplesheet rows require sample, long_reads, short_reads_1, and short_reads_2"
            }

            def meta = [id: row.sample, single_end: false]
            tuple(
                meta,
                file(row.long_reads, checkIfExists: true),
                [
                    file(row.short_reads_1, checkIfExists: true),
                    file(row.short_reads_2, checkIfExists: true)
                ]
            )
        }

    short_reads_ch = samples_ch.map { meta, long_reads, short_reads -> tuple(meta, short_reads, []) }
    long_reads_ch  = samples_ch.map { meta, long_reads, short_reads -> tuple(meta, [], long_reads) }

    FASTP(short_reads_ch, false, false, false)
    FILTLONG(long_reads_ch)

    filtered_reads_ch = FASTP.out.reads
        .join(FILTLONG.out.reads)
        .map { meta, short_reads, long_reads -> tuple(meta, short_reads, long_reads) }

    fastqc_input_ch = FASTP.out.reads
        .mix(FILTLONG.out.reads)
        .map { meta, reads -> tuple(meta, reads) }
    FASTQC(fastqc_input_ch)

    emit:
    reads = filtered_reads_ch
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
}
