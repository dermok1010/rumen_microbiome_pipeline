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

qiime tools export \
  --input-path $BASE_DIR/scripts/ct_table.qza \
  --output-path $BASE_DIR/results/

