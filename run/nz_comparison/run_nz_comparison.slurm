#!/bin/bash
#SBATCH --job-name=nf_nz_cmp
#SBATCH --output=nf_nz_cmp_%j.out
#SBATCH --error=nf_nz_cmp_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#
# Stand-alone 16S run of the NZ_comparison samples. Submit from the run
# directory (NZ_RUN_DIR), NOT from the pipeline checkout, so this run has
# its own .nextflow/ cache and work/ dir and can't hijack the combined
# run's `-resume`:
#   PIPE_DIR=/path/to/rumen_microbiome_pipeline sbatch "$PIPE_DIR/run/nz_comparison/run_nz_comparison.slurm"
set -euo pipefail
: "${PIPE_DIR:?Set PIPE_DIR to the HPC checkout of rumen_microbiome_pipeline}"

module load nextflow
module load apptainer

cd "$SLURM_SUBMIT_DIR"
[ "$(realpath "$SLURM_SUBMIT_DIR")" != "$(realpath "$PIPE_DIR")" ] || {
  echo "ERROR: submit from a separate run directory, not the pipeline checkout." >&2; exit 1; }

PARAMS="$SLURM_SUBMIT_DIR/nz_comparison.params.yaml"
test -f "$PARAMS" || { echo "ERROR: copy and fill nz_comparison.params.yaml into $SLURM_SUBMIT_DIR first." >&2; exit 1; }
grep -v '^ *#' "$PARAMS" | grep -q 'PIPE_DIR\|NZ_RUN_DIR' && { echo "ERROR: placeholders left in $PARAMS." >&2; exit 1; }
CLASSIFIER=$(sed -n 's/^classifier: *"\(.*\)"/\1/p' "$PARAMS")
test -f "$CLASSIFIER" || { echo "ERROR: classifier not found: $CLASSIFIER" >&2; exit 1; }
ls "$PIPE_DIR"/manifest_files/nz_comparison/manifest_*.tsv >/dev/null || {
  echo "ERROR: no manifests -- run make_manifests.py --batch NZcmp first." >&2; exit 1; }

# nextflow.config's hpc profile expects the image at <launch dir>/rumen_microbiome_pipeline.sif.
test -e rumen_microbiome_pipeline.sif || ln -s "$PIPE_DIR/rumen_microbiome_pipeline.sif" .
test -f rumen_microbiome_pipeline.sif || { echo "ERROR: $PIPE_DIR/rumen_microbiome_pipeline.sif missing." >&2; exit 1; }
mkdir -p data/raw   # bound unconditionally by the hpc profile

export NXF_SYNTAX_PARSER=v1
# Record exactly what ran, next to the results.
git -C "$PIPE_DIR" rev-parse HEAD > run_commit.txt
cp "$PARAMS" run_params.yaml

# main.nf's own directory supplies nextflow.config (Nextflow loads the
# project-dir config automatically for a local script path).
nextflow run "$PIPE_DIR/main.nf" \
  -entry full_run -profile hpc -params-file "$PARAMS" -resume
