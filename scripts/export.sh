#!/bin/bash

# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Export the QIIME2 feature table (BIOM format)
qiime tools export \
  --input-path "$BASE_DIR/outputs/ct_table.qza" \
  --output-path "$BASE_DIR/results/"

# Convert BIOM to TSV
biom convert \
  --input-fp "$BASE_DIR/results/feature-table.biom" \
  --output-fp "$BASE_DIR/results/feature-table.tsv" \
  --to-tsv

# Export taxonomy
qiime tools export \
  --input-path "$BASE_DIR/outputs/ct_taxonomy.qza" \
  --output-path "$BASE_DIR/results/"

# Rename taxonomy file for clarity
mv "$BASE_DIR/results/taxonomy.tsv" "$BASE_DIR/results/exported-taxonomy.tsv"

# Completion message
echo " Feature table and taxonomy exported to: $BASE_DIR/results"
