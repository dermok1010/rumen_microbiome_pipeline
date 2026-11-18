#!/bin/bash
# Classify EVERY ASV of the combined 17-flowcell run (1 Oct 2026 table, all host groups) with the
# GTDB r214 515F/806R classifier, for rumen-core's move from SILVA 138 to GTDB r214.
#
# 03_classify_asvs.sh only covered the 505 sheep + 116 Tully beef of the 16S-vs-RE-RRS comparison.
# rumen-core needs every sheep cohort (incl. CT24, CT25, Haggot), Tully 2020 + 2024 and all dairy,
# so this step classifies all ASVs with reads in any non-control sample.
#
# Steps
#  1. Export the representative sequences (outputs/rumen_combined_rep_seqs.qza, 15 Sep run; the
#     1 Oct rerun's rep seqs stayed on the HPC) to FASTA, to list which feature IDs have a sequence.
#  2. 06_full_genus_table.py --write-ids: feature IDs with reads in a non-control sample that also
#     have a sequence -> needed_asvs.tsv; those without a sequence -> asvs_without_rep_seqs.tsv.
#  3. filter-seqs + classify-sklearn (same settings as 03_classify_asvs.sh), export taxonomy.tsv.
#  4. 06_full_genus_table.py builds the genus tables, reconciliation and SILVA->GTDB crosswalk.
#
# Output: outputs/gtdb_r214_full/ (gitignored). Run from the repo root:
#     bash classifier/build/gtdb_r214/05_classify_all_asvs.sh
export OUT=${OUT:-$(cd "$(dirname "$0")/../../.." && pwd)/outputs/gtdb_r214_full}
source "$(dirname "$0")/common.sh"
HERE=$(cd "$(dirname "$0")" && pwd)

if [ ! -s "$OUT/rep_seqs_sep15.fasta" ]; then
    log "export rep seqs to FASTA"
    qiime tools export --input-path "$REPO/outputs/rumen_combined_rep_seqs.qza" --output-path "$OUT/exported_rs"
    mv "$OUT/exported_rs/dna-sequences.fasta" "$OUT/rep_seqs_sep15.fasta" && rmdir "$OUT/exported_rs"
fi
"$PY" "$HERE/06_full_genus_table.py" --write-ids
log "filter rep seqs to $(($(wc -l < "$OUT/needed_asvs.tsv") - 1)) ASVs"
qiime feature-table filter-seqs \
    --i-data "$REPO/outputs/rumen_combined_rep_seqs.qza" \
    --m-metadata-file "$OUT/needed_asvs.tsv" \
    --o-filtered-data "$OUT/all_rep_seqs.qza"
log "classify-sklearn"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-515F-806R-classifier.qza" \
    --i-reads "$OUT/all_rep_seqs.qza" --p-n-jobs "$THREADS" --p-reads-per-batch 2000 \
    --o-classification "$OUT/all_taxonomy_gtdb_r214.qza"
qiime tools export --input-path "$OUT/all_taxonomy_gtdb_r214.qza" --output-path "$OUT/exported_tx"
mv "$OUT/exported_tx/taxonomy.tsv" "$OUT/all_taxonomy_gtdb_r214.tsv" && rmdir "$OUT/exported_tx"
log "taxonomy: $OUT/all_taxonomy_gtdb_r214.tsv"
"$PY" "$HERE/06_full_genus_table.py"
log "done"
