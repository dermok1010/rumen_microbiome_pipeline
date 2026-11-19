# Shared paths and the QIIME2 Docker wrapper for the Greengenes2 2024.09 build. Sourced by the
# step scripts. Mirrors classifier/build/gtdb_r214/common.sh (same image and wrapper).
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
REF=${REF:-$HOME/refs/greengenes2_2024.09}          # downloaded GG2 release files (VM only)
ASSETS=$REPO/classifier/assets/gg2_2024.09          # trained classifier (gitignored)
OUT=${OUT:-$REPO/outputs/gg2_2024.09_full}           # ASV taxonomy + genus tables (gitignored)
GTDB_OUT=$REPO/outputs/gtdb_r214_full                # rep seqs / ASV list from the GTDB r214 run
IMAGE=ghcr.io/dermok1010/rumen-microbiome-pipeline:latest
THREADS=${THREADS:-4}
PY=${PY:-$HOME/envs/rumen_microbiome_pipeline/bin/python}
mkdir -p "$ASSETS/tmp" "$OUT"

qiime() {
    docker run --rm -u "$(id -u):$(id -g)" -e HOME="$ASSETS/tmp" -e TMPDIR="$ASSETS/tmp" \
        -e MPLCONFIGDIR="$ASSETS/tmp" -v "$REPO:$REPO" -v "$REF:$REF:ro" -w "$REPO" "$IMAGE" qiime "$@"
}

log() { echo "[$(date '+%F %T')] $*"; }
