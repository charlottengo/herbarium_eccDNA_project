#!/bin/bash
#SBATCH --job-name=replicon_to_panel
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o /global/home/hpc6076/pacbio_pilot/logs/replicon_to_panel.%j.out
#SBATCH -e /global/home/hpc6076/pacbio_pilot/logs/replicon_to_panel.%j.err
#SBATCH --mail-user=charlotte.ngo@queensu.ca
#SBATCH --mail-type=ALL

set -euo pipefail
log(){ echo "[${SLURM_JOB_ID:-nojob} $(date +'%F %T')] $*"; }

RUN="micromamba run -n flye39"

# ---------- Inputs (override with --export=VAR=...) ----------
REPL=${REPL:-$HOME/pacbio_pilot/data/Amaranthus_palmeri_eccDNA.MT025716.1.fa}
PANEL=${PANEL:-$HOME/pacbio_pilot/work/epsps_scan/circle_panel.fa}
OUTD=${OUTD:-$HOME/pacbio_pilot/work/replicon_to_panel}

PAF="$OUTD/replicon_to_panel.paf"
TSV="$OUTD/replicon_to_panel.filtered.tsv"
BED="$OUTD/replicon_to_panel.bed"

mkdir -p "$OUTD"

log "Replicon: $REPL"
log "Panel:    $PANEL"
log "Out dir:  $OUTD"

[ -s "$REPL" ]  || { log "ERROR: missing replicon FASTA: $REPL"; exit 2; }
[ -s "$PANEL" ] || { log "ERROR: missing circle panel FASTA: $PANEL"; exit 2; }

# ---------- 1) Align replicon to circle panel ----------
log "Running minimap2 (replicon -> circles)..."
$RUN minimap2 -x asm5 -t "${SLURM_CPUS_PER_TASK:-4}" "$REPL" "$PANEL" > "$PAF"

log "minimap2 done. PAF: $PAF"

# Filter good replicon-like hits: L >= 2000 bp, pid >= 0.95
awk '{
  q=$1; qlen=$2; qs=$3; qe=$4; strand=$5;
  t=$6; tlen=$7; ts=$8; te=$9;
  m=$10; L=$11;
  pid = (L>0 ? m/L : 0);
  if (L >= 2000 && pid >= 0.95) {
    printf("%s\t%s\t%d\t%d\t%s\t%.3f\t%d\t%d\t%d\n",
           t,q,ts,te,strand,pid,L,qs,qe);
  }
}' "$PAF" \
  | sort -k1,1 -k3,3n -k6,6nr > "$TSV"

# BED for IGV
awk 'BEGIN{OFS="\t"} {
  circle=$1; q=$2; ts=$3; te=$4; strand=$5; pid=$6; L=$7; qs=$8; qe=$9;
  annot = q"|pid="pid"|len="L"|q:"qs"-"qe"|strand="strand;
  print circle, ts, te, annot;
}' "$TSV" > "$BED"

log "Filtered hits TSV: $TSV"
log "BED for IGV:       $BED"

log "Top hits:"
if [ -s "$TSV" ]; then
  column -t "$TSV" | head
else
  log "No replicon-like hits above thresholds (L>=2000, pid>=0.95)."
fi

log "DONE."
