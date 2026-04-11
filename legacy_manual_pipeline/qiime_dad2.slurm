#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie 
#SBATCH --job-name=qiime_dada2
#SBATCH --output=qiime_dada2_%j.out
#SBATCH --error=qiime_dada2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G                    # Increased memory for paired-end

# Load QIIME2 module
module load qiime2/amp-2024.2

# Activate the conda environment
source activate qiime2-amplicon-2024.2

# Define working directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Run DADA2 denoise-paired
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $BASE_DIR/outputs/paired_end_demux.qza \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 200 \
  --p-n-threads 8 \
  --o-table $BASE_DIR/outputs/table.qza \
  --o-representative-sequences $BASE_DIR/outputs/rep_seqs.qza \
  --o-denoising-stats $BASE_DIR/outputs/denoising_stats.qza \
  --verbose
