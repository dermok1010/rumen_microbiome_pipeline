#!/bin/bash
#SBATCH --job-name=extract_V4
#SBATCH --output=extract_V4_%j.out
#SBATCH --error=extract_V4_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

cd "$SLURM_SUBMIT_DIR"

qiime feature-classifier extract-reads \
  --i-sequences ../assets/silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads ../assets/silva-138-99-515F-806R.qza
