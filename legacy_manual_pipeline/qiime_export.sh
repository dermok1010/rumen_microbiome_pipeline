#!/bin/bash

# Load QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"
EXPORT_DIR="$BASE_DIR/exported"

# Create export directory if it doesn't exist
mkdir -p $EXPORT_DIR

# Export feature table
qiime tools export \
  --input-path $BASE_DIR/scripts/table.qza \
  --output-path $EXPORT_DIR

# Export taxonomy
qiime tools export \
  --input-path $BASE_DIR/results/ct_taxonomy.qza \
  --output-path $EXPORT_DIR

# Convert biom to TSV
biom convert \
  -i $EXPORT_DIR/feature-table.biom \
  -o $EXPORT_DIR/ct_feature-table.tsv \
  --to-tsv
