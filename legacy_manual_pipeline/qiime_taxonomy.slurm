#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_taxonomy
#SBATCH --output=qiime_taxonomy_%j.out
#SBATCH --error=qiime_taxonomy_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G



# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Assign taxonomy using SILVA classifier
qiime feature-classifier classify-sklearn \
  --i-classifier $BASE_DIR/classifier/silva-138-99-341F-806R-classifier.qza \
  --i-reads $BASE_DIR/outputs/rep_seqs.qza \
  --o-classification $BASE_DIR/outputs/ct_taxonomy.qza \
  --p-n-jobs 4

# Summarize taxonomy table for inspection
qiime metadata tabulate \
  --m-input-file $BASE_DIR/outputs/ct_taxonomy.qza \
  --o-visualization $BASE_DIR/outputs/ct_taxonomy.qzv
