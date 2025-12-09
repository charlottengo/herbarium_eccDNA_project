Pilot pipeline for eccDNA and EPSPS-amplicon analysis using published PacBio HiFi data (Carvalho-Moore et al., 2025)

This repository contains a fully documented, modular workflow for detecting and characterizing extrachromosomal circular DNA (eccDNA) and EPSPS amplification structures in Amaranthus palmeri, using published long-read sequencing (PacBio HiFi) as a pilot dataset.

The pipeline is designed to:

Reproduce and validate the eccDNA architecture reported in the literature.

Develop a contig-level “eccDNA panel” that will later be used to detect eccDNA-like signatures in historical herbarium samples (Illumina).

Identify EPSPS locations in the genome and in eccDNA-derived contigs.

Extract depth profiles for EPSPS-bearing contigs to distinguish eccDNA vs. chromosomal amplification.

All steps are implemented using SLURM + micromamba for strict reproducibility.

Repository Structure
herbarium_eccDNA_project/
│
├── bin/                  # All SLURM job scripts used in the workflow
├── docs/                 # Step-by-step modules documenting the pipeline
├── data/                 # Reference EPSPS sequences & replicon (small files only)
├── results/              # Final outputs (depth profiles, EPSPS hits, eccDNA panels)
└── work/                 # Intermediate assembly, mapping, QC, eccDNA detection

Software & Environment

All analyses use isolated micromamba environments to ensure reproducibility.

Install environment
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config append channels defaults
micromamba config set channel_priority strict

micromamba create -n flye39 -y \
    sra-tools seqkit flye minimap2 samtools gzip


For eccDNA detection:

micromamba create -n eccdna -y \
    python=3.10 minimap2 samtools bedtools seqkit pysam

Workflow Overview

Each numbered module in docs/ corresponds to one stage of the pipeline.

01 — Download HiFi Reads & Setup Project

docs/01_download_data.md

Creates a structured project directory.

Downloads the PacBio HiFi dataset (SRR30359588).

Generates basic FASTQ QC with seqkit stats.

Script: bin/sra_fetch.sbatch

02 — Flye Genome Assembly

docs/02_flye_assembly.md

Runs Flye on HiFi reads.

Computes N50, contig counts, per-contig metadata.

Maps HiFi reads back for QC (flagstat, idxstats, coverage).

BUSCO genome completeness.

Scripts:

bin/flye_asm.sbatch

bin/flye_map.sbatch (if separated)

BUSCO summary extraction.

03 — EPSPS & Replicon Mapping

docs/03_epsps_mapping.md

Includes three analytical workflows:

epsps2replicon
Align EPSPS cDNA (FJ861243.1) to the canonical eccDNA replicon (MT025716.1).
→ Confirms replicon structure and EPSPS location.

replicon_to_panel
Align the canonical replicon to the full eccDNA panel.
→ Identifies which circles are “replicon-like”.

scan_EPSPS_genome
Scan the entire Flye assembly for EPSPS hits.
→ Identifies EPSPS-bearing contigs (e.g., contig_631, contig_1406, contig_1407).

Scripts:

bin/epsps2replicon.sbatch

bin/replicon_to_panel.sbatch

bin/scan_EPSPS_genome.sbatch

These analyses revealed the uneven, multi-contig distribution of EPSPS characteristic of eccDNA-derived assemblies.

04 — Depth Profiling

docs/04_depth_profiling.md

Compute per-base coverage across the whole genome.

Extract depth profiles for EPSPS contigs.

Compare depth vs. genome-wide distribution to distinguish:

eccDNA (high, uneven depth)

chromosomal tandem duplication (uniform plateau)

Script:

bin/depth_per_base.sbatch

05 — eccDNA Detection & Circle Panel Construction

docs/05_circle_detection.md

Self-align assembly to identify contigs with head–tail overlaps.

Apply size filters (50 kb–1.5 Mb).

Map HiFi reads to candidates.

Compute coverage fraction per contig.

Validate circularity using junction-spanning reads.

Produce a validated eccDNA panel (FASTA + BED).

Script:

bin/circle_detect.sbatch

helper: bin/circle_qc.py

This procedure recovered 222 high-confidence eccDNA contigs, consistent with large replicon-like circles and additional eccDNA families in the genome.

Scientific Goals of the Pilot

This PacBio re-analysis serves as a controlled test case to:

validate the bioinformatic workflow for eccDNA discovery,

generate a reference eccDNA panel for mapping historical samples,

clarify assembly signatures of eccDNA vs. chromosomal duplications,

learn the depth and structural properties of EPSPS amplifications.

This groundwork enables your PhD Chapters 1 & 2:

determining the timing and origin of EPSPS amplification,

distinguishing eccDNA vs. chromosomal duplication in historical herbarium specimens.

Reproducibility

Every step is:

script-based (in bin/)

environment-controlled (micromamba)

modular (docs/ with numbered steps)

SLURM-compatible for HPC execution

This makes the workflow fully reproducible for committee review and future expansion (e.g., historical Illumina data in Chapter 1).

Future Extensions

Planned modules:

06 — Gene annotation on eccDNA variants

07 — Extracting eccDNA structural variants across populations

08 — Historical Illumina mapping to eccDNA panel

09 — EPSPS amplification trajectories across time

These will be added as the dissertation progresses.
