#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dermot.kelly@teagasc.ie
#SBATCH --job-name=extract_reads
#SBATCH --output=extract_reads_%j.out
#SBATCH --error=extract_reads_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Load and activate QIIME2
module load qiime2/amp-2024.2
source activate qiime2-amplicon-2024.2

# Move into your classifier folder (adjust path if needed)
cd /home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/classifier

# Extract 341F–806R region
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 292 \
  --o-reads silva-138-99-341F-806R.qza
