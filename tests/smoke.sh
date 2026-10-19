#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

printf '@stub_read\nA\n+\nI\n' | gzip -n > "$test_dir/long.fastq.gz"
cp "$test_dir/long.fastq.gz" "$test_dir/r1.fastq.gz"
cp "$test_dir/long.fastq.gz" "$test_dir/r2.fastq.gz"
printf 'sample,long_reads,short_reads_1,short_reads_2\nstub,%s/long.fastq.gz,%s/r1.fastq.gz,%s/r2.fastq.gz\n' \
    "$test_dir" "$test_dir" "$test_dir" > "$test_dir/samplesheet.csv"

NXF_ANSI_LOG=false nextflow run "$project_dir/main.nf" \
    -profile test_stub \
    -stub-run \
    --input "$test_dir/samplesheet.csv" \
    --outdir "$test_dir/results"

test -f "$test_dir/results/polish/stub/stub_final_assembly.fasta"
test -f "$test_dir/results/assembly_qc/stub/stub.mapping_stats.tsv"
grep -q $'Metric\tValue' "$test_dir/results/assembly_qc/stub/stub.mapping_stats.tsv"