# Bacterial Hybrid Assembly Pipeline

A pipeline for bacterial genome long-read assembly, followed by Illumina short read polishing. Primary assembly tool is [Autocycler](https://github.com/rrwick/Autocycler).

## Overview

This pipeline performs end-to-end bacterial genome assembly, from raw reads to a polished assembly with quality metrics.

High-quality nanopore reads are essential -- the pipeline assembles long reads, then polishes with short reads.

## Setup

### Conda install

Navigate to a directory where you wish to clone this repo, then run the following:

```bash
git clone https://github.com/Cmnorris7/bacteria_hybrid_assembly.git
cd ./bacteria_hybrid_assembly/
bash ./setup_environments.sh
```
This setup script will clone the Autocycler repository (a dependency) into the parent directory (i.e., alongside `bacteria_hybrid_assembly`) and create the necessary Conda environments.

## To Run

### Input Requirements

1.  **Create a dedicated directory for your sample's reads.**
2.  Place your raw sequencing files into this directory:
    *   **Nanopore reads:** The filename must contain "nanopore" (e.g., `myisolate_nanopore.fastq.gz`). The file should be gzipped.
        *   The pipeline will automatically determine the **sample name** from the Nanopore filename, using the portion of the name before the first underscore (e.g., `SampleX_nanopore.fastq.gz` will result in `SampleX`).
    *   **Illumina paired-end reads:**
        *   Forward reads filename must contain "R1" (e.g., `myisolate_R1.fastq.gz`).
        *   Reverse reads filename must contain "R2" (e.g., `myisolate_R2.fastq.gz`).
        *   These files should also be gzipped.
    *   These patterns are defined in `config.sh` (`NANOPORE_PATTERN`, `R1_PATTERN`, `R2_PATTERN`).

### Execution

Navigate into the directory containing your Nanopore and Illumina FASTQ files. Then, run the main pipeline script:

```bash
bash /path/to/bacteria_hybrid_assembly/run_bact_assembly.sh
```
(Replace `/path/to/bacteria_hybrid_assembly/` with the actual path to the cloned pipeline repository).

If you run `run_bact_assembly.sh` from a directory whose name does not match the derived `SAMPLE_NAME`, the script will automatically create a subdirectory named after `SAMPLE_NAME`, move the FASTQ files into it, and perform all processing within that subdirectory.

## Pipeline Steps

The main script `run_bact_assembly.sh` orchestrates the following steps:

1.  **Read Quality Control** (`readQC_run.sh`)
    *   Filters Illumina reads with [fastp](https://github.com/OpenGene/fastp):
        *   Filters reads with length < 90bp (configurable via `SHORT_READ_LENGTH` in `config.sh`).
    *   Filters Nanopore reads with [Filtlong](https://github.com/rrwick/Filtlong):
        *   Filters reads with length < 1000bp.
        *   Keeps the best 90% of reads.
        *   Targets a total of 500 Mbp of reads.
    *   Generates [FastQC](https://github.com/s-andrews/FastQC) reports for the filtered reads.
    *   Original raw reads are moved to a `raw_reads/` subdirectory. Filtered reads are placed in `filtered_reads/`.

2.  **Assembly** (`autocycler_run.sh`)
    *   Uses [Autocycler](https://github.com/rrwick/Autocycler) (specifically the `autocycler_full.sh` script from the cloned Autocycler repository) for initial assembly from the filtered long reads.
    *   Autocycler orchestrates multiple assembly tools in parallel to generate an initial consensus assembly, which is output to the `autocycler_out/` directory.

3.  **Assembly Polishing** (`polish_run.sh`)
    *   Reorients the consensus assembly with [dnaapler](https://github.com/gbouras13/dnaapler).
    *   Performs long-read polishing on the reoriented assembly using [Medaka](https://github.com/nanoporetech/medaka) with the model specified in `config.sh` (default: `r1041_e82_400bps_hac_v4.3.0`).
    *   Performs short-read polishing using [Polypolish](https://github.com/rrwick/Polypolish).
    *   Conducts a final round of short-read polishing using [pypolca](https://github.com/gbouras13/pypolca).
    *   A custom Python script (`sort_fasta_add_length.py`) is used to sort contigs by length and append length information to FASTA headers for the final assembly.
    *   Polishing intermediates and the final polished assembly are stored in the `polish/` directory.

4.  **Assembly QC** (`assemblyQC_run.sh`)
    *   Maps the filtered Nanopore reads back to the final polished assembly using [minimap2](https://github.com/lh3/minimap2).
    *   Uses [samtools](https://github.com/samtools) for BAM file processing (sorting, indexing).
    *   Generates detailed coverage and assembly statistics using a custom script (`assemblyQC_coverage_stats.sh`).
    *   Results are stored in the `assembly_qc/` directory, and key statistics are also copied to the main sample directory.

## Output

-   **Final polished assembly:** `polish/${SAMPLE_NAME}_final_assembly.fasta` (where `${SAMPLE_NAME}` is derived from your input file).
-   **Assembly statistics:** `nanopore_map_assembly_stats.tsv` (located in the main sample directory and also in `assembly_qc/`).
-   Detailed logs are stored in the `logs/` directory.

## Customization

Key parameters and environments can be modified in `config.sh`:

-   `MEDAKA_MODEL`: Default is `r1041_e82_400bps_hac_v4.3.0`.
-   `THREADS`: Default is `$(nproc --all)` (uses all available processor cores). You can set a specific number.
-   `SHORT_READ_LENGTH`: Default minimum length for Illumina reads after filtering is 90bp.
-   File patterns for identifying input reads (`NANOPORE_PATTERN`, `R1_PATTERN`, `R2_PATTERN`).
-   Filtering parameters and other tool options can be edited in the respective step scripts (e.g., `readQC_run.sh`, `polish_run.sh`).

## Directory Structure

Within your sample's processing directory, the pipeline will create the following structure:

-   `filtered_reads/`: Quality-filtered Nanopore and Illumina reads.
-   `raw_reads/`: Original input FASTQ files (moved here after initial processing).
-   `autocycler_out/`: Output from the Autocycler assembly step, including the initial consensus assembly.
-   `polish/`: Intermediates from the polishing steps and the final polished assembly.
    -   `polish/tmp/`: Temporary directory used during polishing (removed upon completion).
-   `assembly_qc/`: Output from the assembly quality control step, including BAM files and detailed statistics.
-   `logs/`: Log files for the main pipeline and individual steps.

If the pipeline is run from a directory not matching the sample name, all the above will be within a subdirectory named after the sample (e.g., `SampleX/`).

## Future Work

-   Add test files