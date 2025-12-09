#!/bin/bash
#SBATCH --job-name=depth_per_base
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/depth_per_base.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/depth_per_base.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT="${PROJECT:-$HOME/pacbio_pilot}"
RUN="micromamba run -n flye39"

BAM="${BAM:-$PROJECT/work/mapping_SRR30359588/SRR30359588.hifi.bam}"
OUTD="${OUTD:-$PROJECT/work/depth}"
OUT_TSV="$OUTD/SRR30359588.depth.per_base.tsv.gz"

mkdir -p "$OUTD"

echo "[INFO] BAM  = $BAM"
echo "[INFO] OUTD = $OUTD"

[ -s "$BAM" ] || { echo "[ERROR] Missing BAM: $BAM" >&2; exit 2; }

# Columns: contig, position, depth
echo "[INFO] Running samtools depth..."
$RUN samtools depth -a -@ "${SLURM_CPUS_PER_TASK:-8}" "$BAM" \
  | gzip > "$OUT_TSV"

echo "[INFO] Per-base depth written to: $OUT_TSV"
