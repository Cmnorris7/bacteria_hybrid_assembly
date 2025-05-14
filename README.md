# Bacterial Hybrid Assembly Pipeline

A pipeline for bacterial genome long-read assembly, followed by Illumina short read polishing. Primary assembly tool is [Autocycler](https://github.com/rrwick/Autocycler).

## Overview

This pipeline performs end-to-end bacterial genome assembly, from raw reads to a polished assembly with quality metrics.

High-quality nanopore reads are essential -- pipeline assembles long reads, then polishes with short reads.

## Setup

### Conda install

Navigate to a directory where you wish to clone this repo, then run the following:

```bash
git clone https://github.com/Cmnorris7/bacteria_hybrid_assembly.git
cd .bacteria_hybrid_assembly/
bash ./setup_environments.sh
```

> **_Note_:** Currently configured for use on the NCAH HPC cluster only. All conda environments stored in `/project/shared/miniconda3/envs/`, all singularity containers stored in `/project/scratch/singularity/`

## To Run

### Input

-   A directory containing:
    -   Nanopore reads (must contain "nanopore" in filename)
    -   Illumina paired-end reads (must contain "R1" and "R2" in filenames)

---

From inside the directory containing your reads, run the following command:

```bash
sbatch ${HOME}/git/gitlab/bacteria_hybrid_assembly/run_bact_assembly.sh
```

## Pipeline Steps

1. **Read Quality Control** (`readQC_sbatch.sh`)

    - Filters Illumina reads with [fastp](https://github.com/OpenGene/fastp)
        - Trims adapters
        - Filters reads with quality score < 15
        - Filters reads with length < 140
    - Filters Nanopore reads with [Filtlong](https://github.com/rrwick/Filtlong)
        - Filters reads with length < 1kbp
        - Throws out the worst 10% of reads
        - Removes worst reads until 500 MBP remain
    - Generates [FastQC](https://github.com/s-andrews/FastQC) reports

2. **Assembly** (`autocycler_sbatch.sh`)

    - Uses [Autocycler](https://github.com/rrwick/Autocycler) for initial assembly from long reads

        > **_Note_:** Autocycler spawns ~28 assembly jobs via slurm. Canu takes the longest.

3. **Assembly Polishing** (`polish_sbatch.sh`)

    - Reorients assembly with [dnaapler](https://github.com/gbouras13/dnaapler)
    - Long-read polishing with [Medaka](https://github.com/nanoporetech/medaka) (model: r1041_e82_400bps_hac_v4.3.0)
    - Short-read polishing with [Polypolish](https://github.com/rrwick/Polypolish)
    - Final polishing with [pypolca](https://github.com/gbouras13/pypolca)

4. **Assembly QC** (`assemblyQC_sbatch.sh`)
    - Maps reads back to assembly using [minimap2](https://github.com/lh3/minimap2) and [samtools](https://github.com/samtools)
    - Generates coverage and assembly statistics

## Output

-   Final assembly: `polish/final_assembly.fasta`
-   Assembly statistics: `nanopore_map_assembly_stats.tsv`

## Customization

Key parameters and environments can be modified in `config.sh`:

-   MEDAKA_MODEL: Default is "r1041_e82_400bps_hac_v4.3.0"
-   THREADS: Default is 48
-   File patterns for identifying input reads
-   Filtering parameters and other tool options can be edited in the respective scripts

## Directory Structure

-   `filtered_reads/`: Quality-filtered reads
-   `raw_reads/`: Original input reads
-   `autocycler_out/`: Initial assembly output
-   `polish/`: Polished assembly files
-   `assembly_qc/`: Assembly quality assessment
-   `logs/`: Pipeline execution logs

## Future Work

-   Make script to configure all conda and singularity dependencies on a new system. Add dependency specifics to repo.
