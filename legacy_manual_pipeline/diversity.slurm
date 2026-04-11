#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_diversity
#SBATCH --output=qiime_diversity_%j.out
#SBATCH --error=qiime_diversity_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G



# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"


# Run QIIME2 core-metrics without phylogeny
qiime diversity core-metrics \
  --i-table "$BASE_DIR/scripts/table.qza" \
  --m-metadata-file "$BASE_DIR/scripts/fixed-sample-metadata.tsv" \
  --p-sampling-depth 11000 \
  --p-no-phylogeny \
  --p-n-jobs 8 \
  --output-dir "$BASE_DIR/results/diversity_core_nophylo"
