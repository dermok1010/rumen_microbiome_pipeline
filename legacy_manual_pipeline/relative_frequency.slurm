#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_taxonomy
#SBATCH --output=qiime_taxonomy_%j.out
#SBATCH --error=qiime_taxonomy_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G



# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

 Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Run relative frequency transformation
qiime feature-table relative-frequency \
  --i-table $BASE_DIR/scripts/ct_table.qza \
  --o-relative-frequency-table $BASE_DIR/results/ct_table_relative.qza
