4.1 Data Retrieval & Project Directory Structure

This section describes the directory layout and the SLURM workflow used to download PacBio HiFi reads (SRR30359588) from Carvalho-Moore et al., 2025 for pipeline validation.

4.1.2 Directory Structure

All analyses are organized under a single, clearly named project directory:

PROJECT_DIR="$HOME/pacbio_pilot_CarvalhoMoore2025"
mkdir -p "$PROJECT_DIR"/{bin,data,logs,results,work}

Folder Descriptions
Folder	Description
data/	Raw and processed FASTQ files (SRA → FASTQ)
logs/	SLURM .out and .err logs
bin/	All SLURM scripts and supporting bash pipelines
work/	Intermediate files (assemblies, mappings, coverage outputs)
results/	Summary tables, QC outputs, final files

To reproduce this layout:

mkdir -p "$HOME/pacbio_pilot_CarvalhoMoore2025"/{bin,data,logs,results,work}

4.1.3 Downloading PacBio HiFi Reads (SRR30359588)

A dedicated SLURM script handles download, FASTQ extraction, compression, and QC statistics.

Place this file in:

pacbio_pilot_CarvalhoMoore2025/bin/sra_fetch.sbatch

bin/sra_fetch.sbatch
#!/bin/bash
#SBATCH --job-name=sra_fetch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot_CarvalhoMoore2025/logs/sra_fetch.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot_CarvalhoMoore2025/logs/sra_fetch.%j.err

set -euo pipefail
MAMBA="micromamba run -n flye39"

ACC="SRR30359588"
OUT="$HOME/pacbio_pilot_CarvalhoMoore2025/data"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUT"

echo "Downloading ${ACC} into ${OUT} with ${THREADS} threads"

# 1) Download SRA object
$MAMBA prefetch "$ACC"

# 2) Convert to FASTQ
$MAMBA fasterq-dump "$ACC" \
  -O "$OUT" \
  -t "$OUT" \
  -e "$THREADS"

# 3) Compress FASTQ
gzip "$OUT/${ACC}.fastq"

# 4) Read statistics
$MAMBA seqkit stats -a "$OUT/${ACC}.fastq.gz" \
  | tee "$OUT/${ACC}.seqkit.stats.txt"

echo "Done. FASTQ at: $OUT/${ACC}.fastq.gz"

Submitting the Job

From inside the project directory:

cd "$HOME/pacbio_pilot_CarvalhoMoore2025"
sbatch bin/sra_fetch.sbatch

Outputs

After the job finishes, you will have:

pacbio_pilot_CarvalhoMoore2025/data/
├── SRR30359588.fastq.gz
└── SRR30359588.seqkit.stats.txt


*.fastq.gz → compressed PacBio HiFi reads

*.seqkit.stats.txt → full QC report including read length distribution

