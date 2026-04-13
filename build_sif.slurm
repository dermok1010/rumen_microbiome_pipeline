#!/bin/bash
#SBATCH --job-name=apptainer_pull
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/apptainer_pull.out
#SBATCH --error=logs/apptainer_pull.err

module load apptainer/1.1.9

cd /home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline

export APPTAINER_CACHEDIR=/home/dermot.kelly/Dermot_analysis/apptainer_cache
export APPTAINER_TMPDIR=/home/dermot.kelly/Dermot_analysis/apptainer_cache

mkdir -p $APPTAINER_CACHEDIR

apptainer pull --force rumen_microbiome_pipeline.sif \
docker://ghcr.io/dermok1010/rumen-microbiome-pipeline:latest
