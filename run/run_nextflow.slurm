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

# nextflow.config's hpc profile expects the Apptainer image right here;
# every process fails (and Nextflow retries each one pointlessly, per
# process { errorStrategy = 'retry'; maxRetries = 2 }) without it. Fail in
# seconds instead of burning queue time across every flowcell's retries.
test -f rumen_microbiome_pipeline.sif || {
  echo "ERROR: rumen_microbiome_pipeline.sif not found in $SLURM_SUBMIT_DIR." >&2
  echo "Build/pull it first: APPTAINER_CACHE_ROOT=<scratch dir> sbatch build_sif.slurm" >&2
  exit 1
}

# Nextflow's newer strict syntax parser rejects the classic multiple
# top-level `workflow` blocks + `-entry` pattern this pipeline uses
# ("the -entry option is not supported with the strict parser"). v1
# keeps the classic parser; harmless to set even on older Nextflow
# versions that don't have a strict parser to opt out of.
export NXF_SYNTAX_PARSER=v1

nextflow run main.nf -entry "$ENTRY" -profile hpc -resume
