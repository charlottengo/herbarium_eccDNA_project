#!/usr/bin/env bash
#
# create_flye39_env.sh
# ---------------------
# This script initializes the micromamba environment used for the
# herbarium eccDNA project (Flye assembly, minimap2 mapping, eccDNA detection).
#
# Usage:
#   bash create_flye39_env.sh
#
# This creates an environment named "flye39".
#

set -euo pipefail

echo "=== Configuring micromamba channels ==="
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config append channels defaults
micromamba config set channel_priority strict

echo "=== Creating environment: flye39 ==="
micromamba create -n flye39 -y \
    sra-tools \
    seqkit \
    flye \
    minimap2 \
    samtools \
    gzip

echo "=== Environment created successfully ==="
echo "To use it, run:"
echo '    micromamba run -n flye39 <tool>'
echo ""
echo "You may also define a shortcut alias:"
echo '    export MAMBA="micromamba run -n flye39"'
echo ""
echo "Done."
