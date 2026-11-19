#!/bin/bash
# Train a q2-feature-classifier naive Bayes classifier on the GG2 2024.09 backbone V4 sequences
# (EMP 515F/806R in-silico extract, the same primers as this pipeline), default settings
# (7-mers, 8192 hashed features), no dereplication (class priors as in GG2's own V4 model).
# Taxonomy is truncated at genus (02a_genus_taxonomy.py): the species-level model (28,054
# classes) was killed for memory on the 16 GB VM (exit 137, 2026-11-20). Only genus is used.
# We train our own rather than load GG2's pre-trained 2024.09.backbone.v4.nb.qza because that file
# is a pickle (code) and its scikit-learn version need not match the pipeline image's.
source "$(dirname "$0")/common.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
"$PY" -I "$HERE/02a_genus_taxonomy.py" "$REF" "$ASSETS/backbone_tax_genus.tsv"
qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-path "$ASSETS/backbone_tax_genus.tsv" --output-path "$ASSETS/backbone_tax_genus.qza"
log "fit-classifier-naive-bayes (GG2 2024.09 backbone V4, genus level)"
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$REF/2024.09.backbone.v4.fna.qza" \
    --i-reference-taxonomy "$ASSETS/backbone_tax_genus.qza" \
    --p-classify--chunk-size 2000 \
    --o-classifier "$ASSETS/gg2-2024.09-backbone-v4-genus-nb.qza"
log "classifier: $ASSETS/gg2-2024.09-backbone-v4-genus-nb.qza"
