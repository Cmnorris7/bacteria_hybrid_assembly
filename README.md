# Bacterial Hybrid Assembly

A reproducible Nextflow DSL2 pipeline for bacterial hybrid genome assembly. It filters paired Illumina and Oxford Nanopore reads, builds an ensemble long-read assembly with Autocycler, reorients and polishes the consensus, and reports long-read mapping and coverage statistics.

## Workflow

1. Filter paired short reads with fastp and long reads with Filtlong.
2. generate FastQC reports for filtered reads.
3. Subsample long reads with Autocycler.
4. Assemble every subset with a configurable ensemble.
5. Compress, cluster, trim, resolve, and combine candidates with Autocycler.
6. Reorient with Dnaapler.
7. Polish with Medaka, Polypolish, and PyPolca.
8. Sort contigs and normalize final FASTA headers.
9. Map filtered long reads to the final assembly and generate coverage statistics.

The local profile uses Flye, Miniasm plus Minipolish, and Raven. The SLURM profile additionally enables metaMDBG, Canu, NECAT, and NextDenovo plus NextPolish.

## Requirements

- Linux
- Java 17 or newer
- Nextflow 26.04 or newer
- Docker for `local` runs
- Apptainer and a configured SLURM client for `slurm` runs

All analysis processes run in versioned containers. Local wrappers use immutable image digests; vendored nf-core modules use their recorded Seqera containers.

## Input

Create a CSV with one isolate per row:

```csv
sample,long_reads,short_reads_1,short_reads_2
isolate_01,reads/isolate_01_nanopore.fastq.gz,reads/isolate_01_R1.fastq.gz,reads/isolate_01_R2.fastq.gz
```

Paths may be absolute or relative to the directory where Nextflow is launched. Sample IDs must be unique. See `assets/samplesheet.csv` for a multi-sample example.

## Run Locally

```bash
nextflow run main.nf \
  -profile local \
  --input samplesheet.csv \
  --outdir results \
  --genome_size 5m \
  --medaka_model r1041_e82_400bps_hac_v4.3.0 \
  --long_read_mode nano-hq
```

The local profile is intended for smaller machines and runs three assemblers. Override process resources for a larger workstation in an additional config:

```groovy
// workstation.config
process.resourceLimits = [cpus: 16, memory: '64 GB', time: '48 h']
```

Pass it with `-c workstation.config`.

Resume an interrupted or previously cached run with `-resume`.

## Run On SLURM

```bash
nextflow run main.nf \
  -profile slurm \
  --input samplesheet.csv \
  --outdir results
```

The SLURM profile uses the seven-assembler ensemble and Apptainer. Images should be pre-pulled or available through the site's configured Apptainer cache on offline clusters. Site-specific queue, account, and storage directives can be supplied in an additional config:

```bash
nextflow -C nextflow.config -c site.config run main.nf \
  -profile slurm \
  --input samplesheet.csv
```

## Main Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--input` | required | Input CSV samplesheet |
| `--outdir` | `results` | Published result directory |
| `--genome_size` | `5m` | Expected haploid genome size |
| `--medaka_model` | `r1041_e82_400bps_hac_v4.3.0` | Medaka model matching read chemistry |
| `--long_read_mode` | `nano-hq` | `nano-hq` for R10/high-accuracy reads or `nano-raw` for older/noisier reads |
| `--assemblers` | profile-specific | Comma-separated assembler selection |
| `--autocycler_subsamples` | `4` | Number of long-read subsets |

At least three distinct assemblers are required by the consensus workflow. Full parameter metadata is in `nextflow_schema.json`.

## Results

```text
results/
├── read_qc/<sample>/
│   ├── fastp/
│   ├── fastqc/
│   └── filtlong/
├── assembly/<sample>/
├── polish/<sample>/<sample>_final_assembly.fasta
├── assembly_qc/<sample>/
│   ├── <sample>.bam
│   ├── <sample>.bam.bai
│   └── <sample>.mapping_stats.tsv
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

Nextflow keeps intermediates in `work/` for caching and `-resume`. Remove them only after retaining the desired results and reports.

## Validation

Run the complete containerized stub workflow without biological data:

```bash
bash tests/smoke.sh
```

This checks all 26 process invocations from input QC through final assembly QC. Real-data scientific validation is intentionally separate because assembler performance depends on organism, coverage, and read chemistry.

## Implementation Notes

- nf-core modules are vendored under `modules/nf-core/`.
- Tools unavailable in the module registry have focused wrappers under `modules/local/`.
- The original shell scripts remain under `src/` and at the repository root as legacy reference code. `main.nf` is the supported entry point.
- Raw input files are never modified or moved.
