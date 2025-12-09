# 03 — EPSPS & Replicon Mapping Workflows

This module contains small SLURM workflows to:

1. Validate that the canonical eccDNA replicon contains EPSPS in the expected position.
2. Map the canonical replicon onto your eccDNA circle panel (from the PacBio pilot).
3. Scan the whole Flye assembly for all EPSPS loci (chromosomal + eccDNA-derived).

All scripts live in `bin/` and assume you are using the `flye39` micromamba environment.

---

## 3.1 `bin/epsps2replicon.sbatch`

**Purpose**

Align EPSPS mRNA (FJ861243.1) to the canonical eccDNA replicon (MT025716.1) to confirm:

- EPSPS is present
- its coordinates and orientation
- alignment identity/length

This is a sanity check and a way to generate a small BED file for IGV.

**Script**

```bash
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
Run

bash
Copy code
cd ~/pacbio_pilot
sbatch bin/epsps2replicon.sbatch
Outputs

work/epsps_on_replicon/epsps_to_replicon.paf — raw alignments

.../epsps_to_replicon.filtered.tsv — high-quality EPSPS–replicon hits

.../epsps_on_replicon.bed — EPSPS position on the replicon (for IGV)

3.2 bin/replicon_to_panel.sbatch
Purpose

Align the canonical eccDNA replicon (MT025716.1) to your validated circle panel to see:

which circles are “replicon-like”

what fraction / segments of the replicon they contain

how the replicon is broken up across circles

By default, this assumes you have a circle panel FASTA at:

~/pacbio_pilot/work/epsps_scan/circle_panel.fa

(You can symlink or copy from circle_panel_SRR30359588/out/circle_panel.SRR30359588.fa.)

Script

bash
Copy code
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
Run

bash
Copy code
cd ~/pacbio_pilot
sbatch bin/replicon_to_panel.sbatch
Outputs

work/replicon_to_panel/replicon_to_panel.paf — raw alignments

.../replicon_to_panel.filtered.tsv — which circles carry replicon segments

.../replicon_to_panel.bed — for IGV (replicon coverage across circles)

3.3 bin/scan_EPSPS_genome.sbatch
Purpose

Scan the entire Flye assembly for all EPSPS hits using the EPSPS cDNA as a query, so you can:

identify every contig that carries EPSPS copies

get basic coordinates, alignment lengths, and % identity

generate a BED file for visualization / depth extraction

This is what you used to identify contig_631, contig_1406, contig_1407, etc.

Script

bash
Copy code
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
Run

bash
Copy code
cd ~/pacbio_pilot
sbatch bin/scan_EPSPS_genome.sbatch
Outputs

work/epsps_scan_genome/epsps_to_genome.paf — all EPSPS hits

.../epsps_to_genome.filtered.tsv — cleaned list of EPSPS loci

.../epsps_hits_genome.bed — IGV/JBrowse visualization

