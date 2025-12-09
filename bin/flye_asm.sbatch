#!/bin/bash
#SBATCH --job-name=flye_asm
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/flye_asm.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/flye_asm.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
MAMBA="micromamba run -n flye39"

READS="${READS:-$HOME/pacbio_pilot/data/SRR30359588.fastq.gz}"
OUT="${OUT:-$HOME/pacbio_pilot/work/flye_SRR30359588}"
THREADS="${SLURM_CPUS_PER_TASK:-24}"

echo "READS   = $READS"
echo "OUT     = $OUT"
echo "THREADS = $THREADS"

if [ ! -s "$READS" ]; then
  echo "ERROR: missing FASTQ: $READS"
  exit 1
fi

mkdir -p "$OUT"

echo "Running Flye..."
$MAMBA flye \
  --pacbio-hifi "$READS" \
  -o "$OUT" \
  -t "$THREADS"

echo "Flye assembly completed in $OUT"
