#!/bin/bash
#SBATCH --job-name=circle_detect
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/circle_detect.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/circle_detect.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=BEGIN,END,FAIL

set -Eeuo pipefail

: "${ACC:=SRR30359588}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-16}}"
: "${FA:=$HOME/pacbio_pilot/work/flye_SRR30359588/assembly.fasta}"
: "${READS:=$HOME/pacbio_pilot/data/SRR30359588.fastq.gz}"
: "${ROOT:=$HOME/pacbio_pilot/circle_panel_${ACC}}"

MAPS="${ROOT}/maps"
PAN="${ROOT}/pan"
OUT="${ROOT}/out"
mkdir -p "$MAPS" "$PAN" "$OUT"

export PATH="$HOME/bin:$PATH"
eval "$($HOME/bin/micromamba shell hook -s bash)"
micromamba activate eccdna

echo "[INFO] ACC=$ACC"
echo "[INFO] FA=$FA"

# 1) detect head–tail overlaps
minimap2 -x asm10 -r2k -g5k "$FA" "$FA" > "$MAPS/${ACC}.self.paf"

awk '($1==$6)&&($3<2000)&&($4>$2-2000)&&($8<2000)&&($9>$7-2000){print $1}' \
  "$MAPS/${ACC}.self.paf" | sort -u > "$PAN/${ACC}.circular.ids"

# 2) size filter
seqkit grep -f "$PAN/${ACC}.circular.ids" "$FA" \
  | seqkit seq -m 50000 -M 1500000 > "$PAN/candidate_circles.fa"

# 3) map reads → coverage fraction
minimap2 -ax map-hifi -t "$THREADS" "$PAN/candidate_circles.fa" "$READS" \
  > "$MAPS/${ACC}.circles.sam"

samtools sort -@ "$THREADS" -o "$MAPS/${ACC}.circles.bam" "$MAPS/${ACC}.circles.sam"
samtools index "$MAPS/${ACC}.circles.bam"

bedtools genomecov -ibam "$MAPS/${ACC}.circles.bam" -d \
 | awk '{cov[$1]+=($3>0); len[$1]++} END{for(c in len){printf "%s\t%.3f\n",c,cov[c]/len[c]}}' \
 > "$MAPS/${ACC}.covfrac.txt"

# 4) doubled contigs → junction-spanners
awk '
  /^>/ {h=$0; next}
  {s[h]=s[h]$0}
  END{for(h in s){n=substr(h,2); seq=s[h]; print ">"n"_dbl"; print seq""seq}}
' "$PAN/candidate_circles.fa" > "$PAN/candidate_circles_dbl.fa"

minimap2 -ax map-hifi -t "$THREADS" "$PAN/candidate_circles_dbl.fa" "$READS" \
  > "$MAPS/${ACC}.circles_dbl.sam"

samtools sort -@ "$THREADS" -o "$MAPS/${ACC}.circles_dbl.bam" "$MAPS/${ACC}.circles_dbl.sam"
samtools index "$MAPS/${ACC}.circles_dbl.bam"

# 5) QC + filter
python ~/pacbio_pilot/bin/circle_qc.py
awk 'NR>1 && $3>=0.80 && $4>=2 {print $1}' "$OUT/circle_qc.tsv" > "$OUT/keep.ids"

seqkit grep -f "$OUT/keep.ids" "$PAN/candidate_circles.fa" > "$OUT/circle_panel.${ACC}.fa"
