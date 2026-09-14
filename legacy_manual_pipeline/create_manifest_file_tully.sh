#!/bin/bash

# -----------------------------------------
# Create QIIME2 manifest for Tully_16S data
# Keeps full suffix (e.g. 10105_S72)
# -----------------------------------------

RAW_DIR="/data/BioScience/primary/R2002_methanepredict/NZ_method_comparison/Tully_16S"

OUT_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/manifest_files"
mkdir -p "$OUT_DIR"

MANIFEST="$OUT_DIR/tully_manifest.tsv"

# Header
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

# Loop through forward reads
for f in "$RAW_DIR"/*_R1_001.fastq.gz; do

    base=$(basename "$f")

    # Keep full suffix, remove only lane/read ending
    sample=$(echo "$base" | sed 's/_L001_R1_001.fastq.gz//')

    fwd="$RAW_DIR/$base"
    rev="$RAW_DIR/$(echo "$base" | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')"

    # Only write if reverse pair exists
    if [[ -f "$rev" ]]; then
        echo -e "$sample\t$fwd\t$rev" >> "$MANIFEST"
    else
        echo "Missing pair for $base"
    fi

done

echo "Manifest created: $MANIFEST"
echo "Number of samples:"
tail -n +2 "$MANIFEST" | wc -l
