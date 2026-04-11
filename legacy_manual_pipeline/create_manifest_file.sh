#!/bin/bash

# Path to FASTQ files
RAW_DIR="/home/dermot.kelly/Dermot_primary/microbiome_data/CT_lambs_microbiome_raw_files"

# Manifest output
OUT_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/manifest_files"
mkdir -p "$OUT_DIR"
MANIFEST="$OUT_DIR/ct_manifest.tsv"

# Write header
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

# Loop through _1.fastq.gz files to get sample IDs
for f in "$RAW_DIR"/*_1.fastq.gz; do
  sample=$(basename "$f" _1.fastq.gz)
  fwd="$RAW_DIR/${sample}_1.fastq.gz"
  rev="$RAW_DIR/${sample}_2.fastq.gz"
  echo -e "$sample\t$fwd\t$rev" >> "$MANIFEST"
done

echo "✅ Wide-form manifest created: $MANIFEST"
