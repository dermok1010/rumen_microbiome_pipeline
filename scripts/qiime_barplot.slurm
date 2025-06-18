#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_barplot
#SBATCH --output=qiime_barplot_%j.out
#SBATCH --error=qiime_barplot_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Create taxonomy barplot
qiime taxa barplot \
  --i-table $BASE_DIR/scripts/ct_table.qza \
  --i-taxonomy $BASE_DIR/results/taxonomy.qza \
  --m-metadata-file $BASE_DIR/manifest_files/ct_manifest.tsv \
  --o-visualization $BASE_DIR/results/taxa-bar-plots.qzv
