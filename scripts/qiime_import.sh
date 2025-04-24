#!/bin/bash

# Activate QIIME2
#source activate qiime2-amplicon-2024.2

# Import FASTQ using manifest
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/manifest_files/manifest_goat.csv \
  --output-path single_end_demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize for visual inspection
qiime demux summarize \
  --i-data single_end_demux.qza \
  --o-visualization single_end_demux.qzv
