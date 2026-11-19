#!/bin/bash
# GTDB r214 taxonomy for the V3-V4 ASVs of the 91 phenotyped INZAC Irish hoggets (BSI comparison run,
# Stud=NZIRL sheep), so the hoggets can join the genus-level 16S vs RE-RRS comparison in
# micro_prediction_NZ (re_rrs_comparison/diversity/) as a fourth sheep population. A copy of
# 13_nz_sheep_classify.sh for the hoggets: same 341F/805R classifier, same forward read orientation,
# same rep-seq file (the hoggets' V3-V4 ASVs come from the same DADA2 run as the NZ lambs'). ASVs and
# libraries are the ones in re_rrs_comparison/prep_inputs/hogget_16s_asv/counts.tsv (one library per
# hogget). Needs outputs/region_comparison_v3v4_v4/rep_seqs_V3V4_all.qza from 01_classify.sh.
# Run from the repo root. PIPELINE_DATA_REPO (default: this repo) is the checkout whose classifier
# assets and outputs/ are used, for running from a second worktree of the repo.
source "$(dirname "$0")/../classifier/build/gtdb_r214/common.sh"
REPO=${PIPELINE_DATA_REPO:-$REPO}
ASSETS=$REPO/classifier/assets/gtdb_r214
OUT=$REPO/outputs/region_comparison_v3v4_v4
export PIPELINE_DATA_REPO=$REPO
"$PY" "$(dirname "$0")/16_hogget_genus.py" --write-ids
log "hoggets: filter rep seqs to $(($(wc -l < "$OUT/asvs_INZAC_hoggets.tsv") - 1)) ASVs"
qiime feature-table filter-seqs --i-data "$OUT/rep_seqs_V3V4_all.qza" --m-metadata-file "$OUT/asvs_INZAC_hoggets.tsv" \
    --o-filtered-data "$OUT/rep_seqs_INZAC_hoggets.qza"
log "hoggets: classify-sklearn with gtdb-r214-341F-805R"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-341F-805R-classifier.qza" \
    --i-reads "$OUT/rep_seqs_INZAC_hoggets.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 --p-read-orientation same \
    --o-classification "$OUT/taxonomy_INZAC_hoggets_gtdb_r214.qza"
qiime tools export --input-path "$OUT/taxonomy_INZAC_hoggets_gtdb_r214.qza" --output-path "$OUT/exp_INZAC_hoggets"
mv "$OUT/exp_INZAC_hoggets/taxonomy.tsv" "$OUT/taxonomy_INZAC_hoggets_gtdb_r214.tsv" && rmdir "$OUT/exp_INZAC_hoggets"
"$PY" "$(dirname "$0")/16_hogget_genus.py"
log "done"
