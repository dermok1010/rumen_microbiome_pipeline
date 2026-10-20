#!/bin/bash
#SBATCH --job-name=classifier_V3V4
#SBATCH --output=classifier_V3V4_%j.out
#SBATCH --error=classifier_V3V4_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#
# Only needed if check_reads.py shows the NZ_comparison reads carry the
# 341F/805R (V3-V4) primers. Builds a SILVA 138.99 classifier for that
# amplicon, same recipe as classifier/build/ (no reference truncation --
# the legacy 341F script truncated references to 292 bp, shorter than the
# ~420 bp merged V3-V4 reads). Submit from classifier/assets/:
#   cd classifier/assets && sbatch ../../run/nz_comparison/build_classifier_341F_805R.slurm
set -euo pipefail

module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

cd "$SLURM_SUBMIT_DIR"
test -f silva-138-99-seqs.qza && test -f silva-138-99-tax.qza || {
  echo "ERROR: run from classifier/assets/ with silva-138-99-{seqs,tax}.qza present." >&2; exit 1; }

qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 4 \
  --o-reads silva-138-99-341F-805R.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-341F-805R.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-341F-805R-classifier.qza
