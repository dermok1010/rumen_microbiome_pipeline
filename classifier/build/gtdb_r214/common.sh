# Shared paths and the QIIME2 Docker wrapper for the GTDB r214 build. Sourced by the step scripts.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
ASSETS=$REPO/classifier/assets/gtdb_r214
OUT=$REPO/outputs/gtdb_r214_reclassification
IMAGE=ghcr.io/dermok1010/rumen-microbiome-pipeline:latest
THREADS=${THREADS:-4}
PY=${PY:-$HOME/envs/rumen_microbiome_pipeline/bin/python}  # pandas, for 04_genus_tables.py
mkdir -p "$ASSETS/tmp" "$OUT"

# Run qiime inside the pipeline image as the VM user, with the repo bind-mounted at its own path
# (as the docker_local profile does). TMPDIR on disk, not the container's /tmp.
qiime() {
    docker run --rm -u "$(id -u):$(id -g)" -e HOME="$ASSETS/tmp" -e TMPDIR="$ASSETS/tmp" \
        -e MPLCONFIGDIR="$ASSETS/tmp" -v "$REPO:$REPO" -w "$REPO" "$IMAGE" qiime "$@"
}

log() { echo "[$(date '+%F %T')] $*"; }
