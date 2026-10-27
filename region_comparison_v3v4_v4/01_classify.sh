#!/bin/bash
# GTDB r214 taxonomy for the comparison ASVs of both regions, each with its own primer-matched
# classifier. V3V4_REPSEQS = the rep seqs .qza from Dermot's 2026-08-30 HPC run of BSI's samples.
# Run from the repo root.
# Read orientation is set to "same" (forward): with the default "auto", classify-sklearn guessed
# reverse-complement for the V3-V4 ASVs from a few atypical leading sequences and left ~95% at
# domain level; forward was confirmed from conserved 16S motifs (V3 start at position ~3).
source "$(dirname "$0")/../classifier/build/gtdb_r214/common.sh"
OUT=$REPO/outputs/region_comparison_v3v4_v4
mkdir -p "$OUT"
"$PY" "$(dirname "$0")/02_genus_tables.py" --write-ids

classify() {  # classify <region> <rep seqs .qza> <classifier name>
    local region=$1 seqs=$2 name=$3
    log "$region: filter rep seqs to $(($(wc -l < "$OUT/asvs_$region.tsv") - 1)) ASVs"
    qiime feature-table filter-seqs --i-data "$seqs" --m-metadata-file "$OUT/asvs_$region.tsv" \
        --o-filtered-data "$OUT/rep_seqs_$region.qza"
    log "$region: classify-sklearn with gtdb-r214-$name"
    qiime feature-classifier classify-sklearn \
        --i-classifier "$ASSETS/gtdb-r214-$name-classifier.qza" \
        --i-reads "$OUT/rep_seqs_$region.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 --p-read-orientation same \
        --o-classification "$OUT/taxonomy_${region}_gtdb_r214.qza"
    qiime tools export --input-path "$OUT/taxonomy_${region}_gtdb_r214.qza" --output-path "$OUT/exp_$region"
    mv "$OUT/exp_$region/taxonomy.tsv" "$OUT/taxonomy_${region}_gtdb_r214.tsv" && rmdir "$OUT/exp_$region"
}

[ -s "$OUT/taxonomy_V4_gtdb_r214.tsv" ] || classify V4 "$REPO/outputs/rumen_combined_rep_seqs.qza" 515F-806R
if [ -n "${V3V4_REPSEQS:-}" ]; then
    # the container only mounts the repo, so bring the HPC file inside it first
    cp "$V3V4_REPSEQS" "$OUT/rep_seqs_V3V4_all.qza"
    classify V3V4 "$OUT/rep_seqs_V3V4_all.qza" 341F-805R
else
    log "V3V4_REPSEQS not set: V3-V4 classification skipped"
fi
log "done"
