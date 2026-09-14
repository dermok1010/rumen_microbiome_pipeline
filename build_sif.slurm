#!/bin/bash
#SBATCH --job-name=apptainer_pull
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/apptainer_pull.out
#SBATCH --error=logs/apptainer_pull.err
set -euo pipefail

module load apptainer/1.1.9
cd "$SLURM_SUBMIT_DIR"

: "${APPTAINER_CACHE_ROOT:?Set APPTAINER_CACHE_ROOT to a scratch dir with space for the image cache}"
export APPTAINER_CACHEDIR="$APPTAINER_CACHE_ROOT"
export APPTAINER_TMPDIR="$APPTAINER_CACHE_ROOT"
mkdir -p "$APPTAINER_CACHEDIR"

apptainer pull --force rumen_microbiome_pipeline.sif \
  docker://ghcr.io/dermok1010/rumen-microbiome-pipeline:latest
