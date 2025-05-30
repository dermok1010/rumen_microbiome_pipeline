#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=qiime_summary
#SBATCH --output=qiime_summary%j.out
#SBATCH --error=qiime_summary%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# qiime should already be loaded, uncomment if not
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Set base directory
BASE_DIR="/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline"

# Summarise ASV table
qiime feature-table summarize \
  --i-table $BASE_DIR/scripts/table.qza \
  --o-visualization $BASE_DIR/results/ct_table.qzv

# Summarise representative sequences
qiime feature-table tabulate-seqs \
  --i-data $BASE_DIR/scripts/rep-seqs.qza \
  --o-visualization $BASE_DIR/results/ct_rep-seqs.qzv

# Summarise DADA2 stats
qiime metadata tabulate \
  --m-input-file $BASE_DIR/scripts/denoising_stats.qza \
  --o-visualization $BASE_DIR/results/ct_stats.qzv
