#!/bin/bash
# Extract the 515F/806R region (same primers as ../extract_reads_515F_806R.slurm), dereplicate
# sequences with identical taxonomy, and train the naive-Bayes classifier.
source "$(dirname "$0")/common.sh"
log "extract-reads 515F/806R"
qiime feature-classifier extract-reads \
    --i-sequences "$ASSETS/gtdb-r214-ssu-seqs.qza" \
    --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT \
    --p-n-jobs "$THREADS" \
    --o-reads "$ASSETS/gtdb-r214-515F-806R-seqs.qza"
log "dereplicate (uniq)"
qiime rescript dereplicate \
    --i-sequences "$ASSETS/gtdb-r214-515F-806R-seqs.qza" \
    --i-taxa "$ASSETS/gtdb-r214-ssu-tax.qza" --p-mode uniq --p-threads "$THREADS" \
    --o-dereplicated-sequences "$ASSETS/gtdb-r214-515F-806R-seqs-derep.qza" \
    --o-dereplicated-taxa "$ASSETS/gtdb-r214-515F-806R-tax-derep.qza"
log "fit-classifier-naive-bayes"
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$ASSETS/gtdb-r214-515F-806R-seqs-derep.qza" \
    --i-reference-taxonomy "$ASSETS/gtdb-r214-515F-806R-tax-derep.qza" \
    --o-classifier "$ASSETS/gtdb-r214-515F-806R-classifier.qza"
log "done"
