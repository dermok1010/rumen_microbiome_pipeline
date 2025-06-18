#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=train_classifier
#SBATCH --output=train_classifier_%j.out
#SBATCH --error=train_classifier_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G

# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Navigate to classifier directory
cd /home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/classifier

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-341F-806R.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-341F-806R-classifier.qza \
