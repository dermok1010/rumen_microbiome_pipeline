#!/bin/bash
#SBATCH --job-name=qiime_import
#SBATCH --output=qiime_import_%j.out
#SBATCH --error=qiime_import_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Load QIIME2 environment (adjust if using conda/mamba modules)
source activate qiime2-amplicon-2024.2

# Import FASTQ using manifest
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/manifest_files/ct_manifest.tsv \
  --output-path paired_end_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Summarize for visual inspection
qiime demux summarize \
  --i-data paired_end_demux.qza \
  --o-visualization paired_end_demux.qzv
