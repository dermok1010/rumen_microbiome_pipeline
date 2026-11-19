#!/bin/bash
# Classify every ASV with reads in a non-control sample (the 83,690 ASVs of
# outputs/gtdb_r214_full/all_rep_seqs.qza, built by gtdb_r214/05_classify_all_asvs.sh) with the
# GG2 2024.09 V4 classifier, same classify-sklearn settings as the GTDB r214 run (confidence 0.7
# default), then build genus tables with 04_genus_table.py.
source "$(dirname "$0")/common.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
log "classify-sklearn"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gg2-2024.09-backbone-v4-genus-nb.qza" \
    --i-reads "$GTDB_OUT/all_rep_seqs.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 \
    --o-classification "$OUT/all_taxonomy_gg2_2024.09.qza"
qiime tools export --input-path "$OUT/all_taxonomy_gg2_2024.09.qza" --output-path "$OUT/exported_tx"
mv "$OUT/exported_tx/taxonomy.tsv" "$OUT/all_taxonomy_gg2_2024.09.tsv" && rmdir "$OUT/exported_tx"
"$PY" "$HERE/04_genus_table.py"
log "done"
