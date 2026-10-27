#!/bin/bash
# V4-core comparison: trim both regions to the same 251 bp V4 core (04_trim_to_v4_core.py),
# classify every core once with the GTDB r214 V4 classifier, build genus tables, and rerun the
# comparison on them (03_compare.R with CORE=1, outputs in outputs/region_comparison_v3v4_v4/core/).
# Run from the repo root after 01_classify.sh.
source "$(dirname "$0")/../classifier/build/gtdb_r214/common.sh"
OUT=$REPO/outputs/region_comparison_v3v4_v4
"$PY" "$(dirname "$0")/04_trim_to_v4_core.py"
log "import and classify $(grep -c '>' "$OUT/core_seqs.fasta") core sequences"
qiime tools import --type 'FeatureData[Sequence]' --input-path "$OUT/core_seqs.fasta" \
    --output-path "$OUT/core_seqs.qza"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-515F-806R-classifier.qza" --i-reads "$OUT/core_seqs.qza" \
    --p-n-jobs "$THREADS" --p-reads-per-batch 2000 --p-read-orientation same \
    --o-classification "$OUT/core_taxonomy_gtdb_r214.qza"
qiime tools export --input-path "$OUT/core_taxonomy_gtdb_r214.qza" --output-path "$OUT/exp_core"
mv "$OUT/exp_core/taxonomy.tsv" "$OUT/core_taxonomy_gtdb_r214.tsv" && rmdir "$OUT/exp_core"
"$PY" "$(dirname "$0")/04_trim_to_v4_core.py" --genus
CORE=1 micromamba run -p "$HOME/envs/micro_prediction" Rscript "$(dirname "$0")/03_compare.R"
log "done"
