#!/bin/bash
#SBATCH --job-name=scan_EPSPS_genome
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/scan_epsps_genome.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/scan_epsps_genome.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=ALL

set -euo pipefail

RUN="micromamba run -n flye39"

REF=${REF:-$HOME/pacbio_pilot/data/EPSPS.FJ861243.1.fa}
ASM=${ASM:-$HOME/pacbio_pilot/work/flye_SRR30359588/assembly.fasta}
OUTD=${OUTD:-$HOME/pacbio_pilot/work/epsps_scan_genome}

PAF="$OUTD/epsps_to_genome.paf"
TSV="$OUTD/epsps_to_genome.filtered.tsv"
BED="$OUTD/epsps_hits_genome.bed"

mkdir -p "$OUTD"

echo "[INFO] REF  = $REF"
echo "[INFO] ASM  = $ASM"
echo "[INFO] OUTD = $OUTD"

[ -s "$REF" ] || { echo "[ERROR] Missing EPSPS REF: $REF" >&2; exit 2; }
[ -s "$ASM" ] || { echo "[ERROR] Missing assembly: $ASM" >&2; exit 2; }

echo "[INFO] Running minimap2 (EPSPS -> genome)..."
$RUN minimap2 -x splice:hq -t "${SLURM_CPUS_PER_TASK:-4}" "$REF" "$ASM" > "$PAF"

echo "[INFO] Filtering hits..."
# PAF: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq ...
awk '{
  q=$1;  qlen=$2;  qs=$3;  qe=$4;
  t=$6;  ts=$8;   te=$9;
  m=$10; L=$11;
  pid = (L>0 ? m/L : 0);
  if (L >= 250 && pid >= 0.90) {
    printf("%s\t%d\t%d\t%d\t%s\t%.3f\t%d\t%d\t%d\n",
           q,qlen,qs,qe,t,pid,L,ts,te);
  }
}' "$PAF" \
  | sort -k5,5 -k6,6nr -k8,8n > "$TSV"

# BED for IGV: tname, tstart, tend, label
awk 'BEGIN{OFS="\t"} {print $5,$8,$9,$1" pid="$6" len="$7}' "$TSV" > "$BED"

echo "[INFO] Done."
echo "[INFO] Filtered hits table : $TSV"
echo "[INFO] BED for IGV         : $BED"

echo "[INFO] Top hits (if any):"
head "$TSV" || echo "[INFO] No hits passing filters."
