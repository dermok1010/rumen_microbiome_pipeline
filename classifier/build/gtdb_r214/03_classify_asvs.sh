#!/bin/bash
# Classify the ASVs present in the comparison animals with the GTDB r214 classifier.
# needed_asvs.tsv (written by 04_genus_tables.py --write-ids) lists their feature IDs.
source "$(dirname "$0")/common.sh"
[ -s "$OUT/needed_asvs.tsv" ] || "$PY" "$(dirname "$0")/04_genus_tables.py" --write-ids
log "filter rep seqs to $(($(wc -l < "$OUT/needed_asvs.tsv") - 1)) ASVs"
qiime feature-table filter-seqs \
    --i-data "$REPO/outputs/rumen_combined_rep_seqs.qza" \
    --m-metadata-file "$OUT/needed_asvs.tsv" \
    --o-filtered-data "$OUT/comparison_rep_seqs.qza"
log "classify-sklearn"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-515F-806R-classifier.qza" \
    --i-reads "$OUT/comparison_rep_seqs.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 \
    --o-classification "$OUT/comparison_taxonomy_gtdb_r214.qza"
qiime tools export --input-path "$OUT/comparison_taxonomy_gtdb_r214.qza" --output-path "$OUT/exported"
mv "$OUT/exported/taxonomy.tsv" "$OUT/comparison_taxonomy_gtdb_r214.tsv" && rmdir "$OUT/exported"
log "done: $OUT/comparison_taxonomy_gtdb_r214.tsv"
