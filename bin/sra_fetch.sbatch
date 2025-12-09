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
