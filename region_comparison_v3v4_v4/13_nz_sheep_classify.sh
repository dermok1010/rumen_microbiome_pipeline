#!/bin/bash
# GTDB r214 taxonomy for the V3-V4 ASVs of the 92 NZ Flock 2638 lambs (BSI comparison run), so the
# NZ sheep can join the genus-level diversity comparison with RE-RRS in micro_prediction_NZ
# (re_rrs_comparison/diversity/genus_diversity.R). Same 341F/805R classifier and forward read
# orientation as the dairy V3-V4 step (01_classify.sh). ASVs and libraries are the ones in
# re_rrs_comparison/prep_inputs/nz_16s_asv/counts.tsv (one library per lamb).
# Needs outputs/region_comparison_v3v4_v4/rep_seqs_V3V4_all.qza from 01_classify.sh.
# Run from the repo root.
source "$(dirname "$0")/../classifier/build/gtdb_r214/common.sh"
OUT=$REPO/outputs/region_comparison_v3v4_v4
NZ=/home/dermodkkelly/microbiome_prediction/re_rrs_comparison/prep_inputs/nz_16s_asv/counts.tsv
"$PY" "$(dirname "$0")/14_nz_sheep_genus.py" --write-ids
log "NZ sheep: filter rep seqs to $(($(wc -l < "$OUT/asvs_NZ_sheep.tsv") - 1)) ASVs"
qiime feature-table filter-seqs --i-data "$OUT/rep_seqs_V3V4_all.qza" --m-metadata-file "$OUT/asvs_NZ_sheep.tsv" \
    --o-filtered-data "$OUT/rep_seqs_NZ_sheep.qza"
log "NZ sheep: classify-sklearn with gtdb-r214-341F-805R"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-341F-805R-classifier.qza" \
    --i-reads "$OUT/rep_seqs_NZ_sheep.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 --p-read-orientation same \
    --o-classification "$OUT/taxonomy_NZ_sheep_gtdb_r214.qza"
qiime tools export --input-path "$OUT/taxonomy_NZ_sheep_gtdb_r214.qza" --output-path "$OUT/exp_NZ_sheep"
mv "$OUT/exp_NZ_sheep/taxonomy.tsv" "$OUT/taxonomy_NZ_sheep_gtdb_r214.tsv" && rmdir "$OUT/exp_NZ_sheep"
"$PY" "$(dirname "$0")/14_nz_sheep_genus.py"
log "done"
