#!/bin/bash
#SBATCH --job-name=make_manifests
#SBATCH --output=make_manifests_%j.out
#SBATCH --error=make_manifests_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
set -euo pipefail

# same conda activation you use elsewhere
source /install/software/2024/minpy3/python3.10/etc/profile.d/conda.sh
conda activate /data/Genetics/analysis/R1681_OviSeq/Dermot/envs/dermo_python_env

# run from wherever you sbatch this; writes manifests into ./manifest_files_combined
python3 make_manifests.py "$PWD/manifest_files/combined_run"
