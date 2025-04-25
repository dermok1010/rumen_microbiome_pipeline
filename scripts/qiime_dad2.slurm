#!/bin/bash
#SBATCH --job-name=qiime_dada2
#SBATCH --output=qiime_dada2_%j.out
#SBATCH --error=qiime_dada2_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load the module qiime2
module load qiime2/amp-2024.2

# Activate QIIME2 environment
source activate qiime2-amplicon-2024.2

# Set base working directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Run DADA2 denoising
qiime dada2 denoise-single \
  --i-demultiplexed-seqs $BASE_DIR/scripts/single_end_demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --o-table $BASE_DIR/table.qza \
  --o-representative-sequences $BASE_DIR/rep-seqs.qza \
  --o-denoising-stats $BASE_DIR/stats.qza
