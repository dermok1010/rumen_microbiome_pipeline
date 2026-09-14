#!/bin/bash
#SBATCH --job-name=train_V4
#SBATCH --output=train_V4_%j.out
#SBATCH --error=train_V4_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G

module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

cd "$SLURM_SUBMIT_DIR"

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ../assets/silva-138-99-515F-806R.qza \
  --i-reference-taxonomy ../assets/silva-138-99-tax.qza \
  --o-classifier ../assets/silva-138-99-515F-806R-classifier.qza
