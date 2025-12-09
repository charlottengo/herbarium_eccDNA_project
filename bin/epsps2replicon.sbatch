#!/bin/bash
#SBATCH --job-name=epsps2replicon
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/epsps2replicon.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/epsps2replicon.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=ALL

set -euo pipefail

RUN="micromamba run -n flye39"

REF=${REF:-$HOME/pacbio_pilot/data/EPSPS.FJ861243.1.fa}
REPL=${REPL:-$HOME/pacbio_pilot/data/Amaranthus_palmeri_eccDNA.MT025716.1.fa}
OUTD=${OUTD:-$HOME/pacbio_pilot/work/epsps_on_replicon}

PAF="$OUTD/epsps_to_replicon.paf"
TSV="$OUTD/epsps_to_replicon.filtered.tsv"
BED="$OUTD/epsps_on_replicon.bed"

mkdir -p "$OUTD"

[ -s "$REF" ]  || { echo "Missing EPSPS REF: $REF";  exit 2; }
[ -s "$REPL" ] || { echo "Missing REPLICON: $REPL"; exit 2; }

echo "[INFO] REF  = $REF"
echo "[INFO] REPL = $REPL"

# map cDNA (EPSPS) to eccDNA replicon
$RUN minimap2 -x splice -t "${SLURM_CPUS_PER_TASK:-4}" "$REPL" "$REF" > "$PAF"

# Filter: keep long, high-identity hits (>= 1000 bp and >= 95% identity)
awk '{
  # PAF: tname tlen tstart tend strand qname qlen qstart qend matches alnlen mapq ...
  t=$1; tlen=$2; ts=$3; te=$4; strand=$5;
  q=$6; qlen=$7; qs=$8; qe=$9; matches=$10; alnlen=$11;
  pid = (alnlen>0 ? matches/alnlen : 0);
  if (alnlen >= 1000 && pid >= 0.95) {
    printf("%s\t%d\t%d\t%d\t%s\t%.3f\t%d\t%d\t%d\t%d\n",
           t,tlen,ts,te,strand,pid,alnlen,qlen,qs,qe);
  }
}' "$PAF" \
  | sort -k1,1 -k3,3n -k6,6nr > "$TSV"

# BED-style output for IGV
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$5" pid="$6" len="$7}' "$TSV" > "$BED"

echo "[INFO] Filtered hits (contig tlen tstart tend strand pid alnlen qlen qs qe):"
head "$TSV" || true
echo "[INFO] BED written: $BED"
echo "[INFO] PAF raw:     $PAF"
