#!/bin/bash
#SBATCH --job-name=nf_rumen
#SBATCH --output=nf_%j.out
#SBATCH --error=nf_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#
# Usage: sbatch run_nextflow.slurm <entry>   (entry: import_only | full_run)
set -euo pipefail
ENTRY="${1:?Usage: sbatch run_nextflow.slurm <import_only|full_run>}"

module load nextflow
module load apptainer

cd "$SLURM_SUBMIT_DIR"

nextflow run main.nf -entry "$ENTRY" -profile hpc -resume
