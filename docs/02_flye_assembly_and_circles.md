# Step 2 — Flye Genome Assembly & eccDNA Detection  
_PacBio HiFi pilot (Carvalho-Moore et al., 2025)_

This module documents the end-to-end pipeline for assembling the PacBio HiFi dataset  
(SRR30359588) and detecting high-confidence eccDNA circles.  
All scripts referenced here live in `bin/`.

---

# 2.1 Flye Genome Assembly

## 2.1.1 SLURM script — `bin/flye_asm.sbatch`

```bash
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
Run it:
bash
Copy code
cd ~/pacbio_pilot
sbatch bin/flye_asm.sbatch
**2.1.2 Flye Outputs**
Location: work/flye_SRR30359588/

File	Description
assembly.fasta	primary Flye assembly
assembly_info.txt	per-contig metadata (length, cov, circularity)
params.json	Flye parameters for reproducibility

**2.2 Assembly Quality Checks**
**2.2.1 Inspect contig metadata**
bash
Copy code
column -t work/flye_SRR30359588/assembly_info.txt | less
Typical features observed:

large contigs (~5–7 Mb) with ~14–15× depth

smaller circular contigs (~27–30×) likely organelles or repeats

**2.2.2 Assembly statistics (contigs, N50, total length)**
bash
Copy code
MAMBA="micromamba run -n flye39"
ASM="$HOME/pacbio_pilot/work/flye_SRR30359588/assembly.fasta"

$MAMBA samtools faidx "$ASM"

awk '
BEGIN { sum = 0 }
     { len = $2; sum += len; L[NR] = len }
END {
  n = asort(L)
  half = sum / 2
  acc = 0
  for (i = n; i >= 1; i--) {
      acc += L[i]
      if (acc >= half) { n50 = L[i]; break }
  }
  printf("contigs = %d\n", n)
  printf("total   = %d bp\n", sum)
  printf("N50     = %d bp\n", n50)
}' "$ASM.fai" \
  > work/flye_SRR30359588/assembly.summary.txt
2.3 Mapping HiFi Reads Back to Assembly
Create mapping directory:

bash
Copy code
MAPDIR="$HOME/pacbio_pilot/work/mapping_SRR30359588"
mkdir -p "$MAPDIR"
**2.3.1 Map reads**
bash
Copy code
zcat "$READS" \
  | minimap2 -ax map-hifi -t "$THREADS" "$ASM" - \
  > "$MAPDIR/SRR30359588.sam"
**2.3.2 Sort & index**
bash
Copy code
samtools view -@ "$THREADS" -b "$MAPDIR/SRR30359588.sam" \
  > "$MAPDIR/SRR30359588.unsorted.bam"

samtools sort -@ "$THREADS" \
  -o "$MAPDIR/SRR30359588.hifi.bam" \
  "$MAPDIR/SRR30359588.unsorted.bam"

samtools index "$MAPDIR/SRR30359588.hifi.bam"
**2.3.3 Coverage**
bash
Copy code
samtools coverage "$MAPDIR/SRR30359588.hifi.bam" \
  > "$MAPDIR/SRR30359588.coverage.tsv"
**2.4 Final QC metrics**
Percent mapped
bash
Copy code
samtools flagstat "$MAPDIR/SRR30359588.hifi.bam"
Mean depth
bash
Copy code
awk 'NR>1 {L+=($3-$2+1); C+=($3-$2+1)*$7} END{print C/L}' \
  "$MAPDIR/SRR30359588.coverage.tsv"
Median contig depth
bash
Copy code
awk 'NR>1{print $7}' "$MAPDIR/SRR30359588.coverage.tsv" \
  | sort -n \
  | awk '{a[NR]=$1} END{print (NR%2?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2)}'
**2.5 BUSCO genome completeness**
bash
Copy code
micromamba run -n flye39 busco \
  -i $HOME/pacbio_pilot/work/flye_SRR30359588/assembly.fasta \
  -m genome \
  -l eudicots_odb10 \
  -o SRR30359588_busco \
  -c 24 \
  --out_path $HOME/pacbio_pilot/results/busco
2.6 eccDNA Circle Detection
eccDNA candidates are detected by:

head–tail overlap detection

size filter (50 kb–1.5 Mb)

read mapping & coverage fraction

junction-spanning validation

strict filtering (cov_frac ≥ 0.80, spanners ≥ 2)

**2.6.1 Environment**
bash
Copy code
micromamba create -n eccdna -y \
  python=3.10 minimap2 samtools bedtools seqkit pysam
**2.6.2 SLURM script – bin/circle_detect.sbatch**
bash
Copy code
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
