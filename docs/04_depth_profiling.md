# 04 — Depth Profiling of EPSPS & eccDNA Contigs

This module describes how to:

1. Generate **per-base depth** across the Flye assembly.
2. Summarise depth per contig.
3. Extract depth profiles for **EPSPS-positive contigs** (e.g., contig_631, contig_1406, contig_1407) to compare them to the genome-wide background.

Inputs (from previous modules):

- Assembly: `work/flye_SRR30359588/assembly.fasta`
- Mapped HiFi reads: `work/mapping_SRR30359588/SRR30359588.hifi.bam`
- Coverage summary: `work/mapping_SRR30359588/SRR30359588.coverage.tsv`
- EPSPS hits table: `work/epsps_scan_genome/epsps_to_genome.filtered.tsv`

All commands assume the project lives at:

```bash
PROJECT=~/pacbio_pilot

4.1 Per-base depth for the whole assembly

We use samtools depth to compute per-base coverage across all contigs.

SLURM script: bin/depth_per_base.sbatch
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

Output:

work/depth/SRR30359588.depth.per_base.tsv.gz

columns: contig position depth

**4.2 Contig-level mean depth **

I already computed contig-level mean depth with samtools coverage:

File: work/mapping_SRR30359588/SRR30359588.coverage.tsv

If I need to recompute:

PROJECT=~/pacbio_pilot
BAM="$PROJECT/work/mapping_SRR30359588/SRR30359588.hifi.bam"
OUT="$PROJECT/work/mapping_SRR30359588/SRR30359588.coverage.tsv"

micromamba run -n flye39 samtools coverage "$BAM" > "$OUT"


This file has one row per contig, with columns including:

#rname – contig name

length – contig length

meandepth – mean depth across that contig

I can quickly pull the genome-wide mean and median again:

COV="$PROJECT/work/mapping_SRR30359588/SRR30359588.coverage.tsv"

# Mean of per-base depth (already done earlier, here for reference)
awk 'NR>1 {L+=($3-$2+1); C+=($3-$2+1)*$7} END{print "mean_depth = "C/L}' "$COV"

# Median of contig mean depths
awk 'NR>1{print $7}' "$COV" \
  | sort -n \
  | awk '{a[NR]=$1} END{print "median_contig_depth = " (NR%2?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2)}'

**4.3 Depth profiles for EPSPS-positive contigs**

Now we want to compare EPSPS-bearing contigs such as contig_631, contig_1406, contig_1407 against the genomic background.

4.3.1 Identify EPSPS contigs

From Module 03 (scan_EPSPS_genome.sbatch), you have:

EPS_TSV=~/pacbio_pilot/work/epsps_scan_genome/epsps_to_genome.filtered.tsv
head "$EPS_TSV"


The columns (from the awk in that script) are:

qname (EPSPS)

qlen

qstart

qend

tname (contig)

%id

aln_len

tstart

tend

I can list unique EPSPS-containing contigs:

awk '{print $5}' "$EPS_TSV" | sort -u

**4.3.2 Extract per-base depth for EPSPS contigs
**
Using the per-base depth file:

PROJECT=~/pacbio_pilot
DEP_GZ="$PROJECT/work/depth/SRR30359588.depth.per_base.tsv.gz"
OUTD="$PROJECT/results/depth_profiles"
mkdir -p "$OUTD"

# Example: extract depth for contig_631
zgrep -E '^contig_631\t' "$DEP_GZ" > "$OUTD/contig_631.depth.tsv"

# Similarly for contig_1406, contig_1407
zgrep -E '^contig_1406\t' "$DEP_GZ" > "$OUTD/contig_1406.depth.tsv"
zgrep -E '^contig_1407\t' "$DEP_GZ" > "$OUTD/contig_1407.depth.tsv"


Each file has:

contig position depth

I can plot these in R to show the crazy spike for contig_631 vs. modest depth for 1406/1407.

**4.4 Tiny R snippet for quick depth visualization**

Example: compare mean depth of a few contigs vs genome-wide:

library(dplyr)
library(ggplot2)
library(readr)

project <- "~/pacbio_pilot"
cov <- file.path(project, "work/mapping_SRR30359588/SRR30359588.coverage.tsv")

cov_df <- read_tsv(cov, comment = "#")  # samtools coverage adds a header

# Focus on EPSPS contigs of interest
epsps_contigs <- c("contig_631", "contig_1406", "contig_1407")

epsps_depth <- cov_df %>%
  filter(`#rname` %in% epsps_contigs) %>%
  select(contig = `#rname`, length = length, meandepth = `meandepth`)

# Genome-wide distribution
ggplot(cov_df, aes(x = meandepth)) +
  geom_histogram(bins = 60, fill = "grey80", colour = "grey40") +
  geom_vline(data = epsps_depth,
             aes(xintercept = meandepth, colour = contig),
             linewidth = 1.1) +
  scale_colour_manual(values = c("contig_631" = "purple4",
                                 "contig_1406" = "dodgerblue3",
                                 "contig_1407" = "darkorange3")) +
  labs(x = "Mean depth per contig",
       y = "Number of contigs",
       colour = "EPSPS contig",
       title = "EPSPS-bearing contigs vs. genome-wide depth") +
  theme_bw(base_size = 14)


I can also make a tiny log-scale barplot:

ggplot(epsps_depth,
       aes(x = contig, y = meandepth, fill = contig)) +
  geom_col() +
  scale_y_log10() +
  labs(y = "Mean depth (log10 scale)",
       x = NULL,
       title = "EPSPS contigs show extreme vs. modest amplification") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

**4.5 Interpretation link to eccDNA vs chromosomal**

Once I have these depth profiles, I can interpret them in the context of your thesis:

_eccDNA-like pattern_

contig_631: mean depth >> genome-wide mean (e.g., ~268× vs ~18×)

likely an amplified eccDNA-derived fragment

weakly amplified fragments

contig_1406/1407: modestly above background (~10–11×)

may represent low-copy insertions or partial replicon fragments

_chromosomal duplication_ would look different

a single contig containing EPSPS in a long, uniform-depth plateau, with smooth depth across flanking genomic sequence

This depth module gives me all the raw pieces to show those contrasts clearly.


---

## 📜 SLURM script file to add

Create this as `bin/depth_per_base.sbatch` in your repo:

```bash
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

echo "[INFO] Running samtools depth..."
$RUN samtools depth -a -@ "${SLURM_CPUS_PER_TASK:-8}" "$BAM" \
  | gzip > "$OUT_TSV"

echo "[INFO] Per-base depth written to: $OUT_TSV"
