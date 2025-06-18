#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_collapse
#SBATCH --output=qiime_collapse_%j.out
#SBATCH --error=qiime_collapse_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G



# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set home directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Define paths (adjust if needed)
INPUT_TABLE="$BASE_DIR/results/ct_table_relative.qza"
TAXONOMY="$BASE_DIR/results/ct_taxonomy.qza"
COLLAPSED_TABLE="$BASE_DIR/results/genus_table_relative.qza"
EXPORT_DIR="$BASE_DIR/results/genus_table_tsv"

# Define file paths
RAW_TABLE="$BASE_DIR/scripts/ct_table.qza"
TAXONOMY="$BASE_DIR/results/ct_taxonomy.qza"
COLLAPSED_TABLE="$BASE_DIR/results/genus_table.qza"
RELATIVE_TABLE="$BASE_DIR/results/genus_table_relative.qza"
EXPORT_DIR="$BASE_DIR/results/genus_table_tsv"

# Step 1: Collapse raw table to genus level
qiime taxa collapse \
  --i-table "$RAW_TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --p-level 6 \
  --o-collapsed-table "$COLLAPSED_TABLE" || { echo "Collapse failed"; exit 1; }

# Step 2: Convert to relative frequency
qiime feature-table relative-frequency \
  --i-table "$COLLAPSED_TABLE" \
  --o-relative-frequency-table "$RELATIVE_TABLE" || { echo "RelFreq failed"; exit 1; }

# Step 3: Export to TSV
mkdir -p "$EXPORT_DIR"

qiime tools export \
  --input-path "$RELATIVE_TABLE" \
  --output-path "$EXPORT_DIR" || { echo "Export failed"; exit 1; }

biom convert \
  --input-fp "$EXPORT_DIR/feature-table.biom" \
  --output-fp "$EXPORT_DIR/genus_relative_abundance.tsv" \
  --to-tsv || { echo "BIOM to TSV conversion failed"; exit 1; }
