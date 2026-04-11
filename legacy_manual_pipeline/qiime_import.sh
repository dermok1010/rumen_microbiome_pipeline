#!/bin/bash

# Activate QIIME2
source activate qiime2-amplicon-2024.2

# Define path to manifest file
MANIFEST="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/manifest_files/ct_manifest.csv"

# Define output filenames
OUTPUT_QZA="ct_paired_end_demux.qza"
OUTPUT_QZV="ct_paired_end_demux.qzv"

# Import paired-end FASTQ data using the manifest
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "$OUTPUT_QZA" \
  --input-format PairedEndFastqManifestPhred33V2

# Summarize demultiplexed data
qiime demux summarize \
  --i-data "$OUTPUT_QZA" \
  --o-visualization "$OUTPUT_QZV"

echo "QIIME2 import and summary complete:"
echo "  Imported data artifact: $OUTPUT_QZA"
echo "  Visualization file:     $OUTPUT_QZV"
