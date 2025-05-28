#!/bin/bash

# Path to FASTQ files
RAW_DIR="/home/dermot.kelly/Dermot_primary/microbiome_data"

# Output file
MANIFEST="manifest.csv"

# Header
echo "sample-id,absolute-filepath,direction" > "$MANIFEST"

# Loop through all _1 files and build paired rows
for f in "$RAW_DIR"/*_1.fastq.gz; do
  sample=$(basename "$f" _1.fastq.gz)
  echo "$sample,$RAW_DIR/${sample}_1.fastq.gz,forward" >> "$MANIFEST"
  echo "$sample,$RAW_DIR/${sample}_2.fastq.gz,reverse" >> "$MANIFEST"
done

echo "Manifest created: $MANIFEST"
