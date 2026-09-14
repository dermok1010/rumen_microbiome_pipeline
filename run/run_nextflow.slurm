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

# Nextflow's newer strict syntax parser rejects the classic multiple
# top-level `workflow` blocks + `-entry` pattern this pipeline uses
# ("the -entry option is not supported with the strict parser"). v1
# keeps the classic parser; harmless to set even on older Nextflow
# versions that don't have a strict parser to opt out of.
export NXF_SYNTAX_PARSER=v1

nextflow run main.nf -entry "$ENTRY" -profile hpc -resume
