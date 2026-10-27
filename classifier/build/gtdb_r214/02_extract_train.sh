#!/bin/bash
# Extract the amplicon region, dereplicate sequences with identical taxonomy, and train the
# naive-Bayes classifier.
#   REGION=V4   (default) 515F/806R, the same primers as ../extract_reads_515F_806R.slurm
#   REGION=V3V4 341F/805R, for BSI's NZ_comparison run (DADA2 trim_left 17/21 = these primer lengths)
source "$(dirname "$0")/common.sh"
case "${REGION:-V4}" in
    V4)   NAME=515F-806R; F=GTGCCAGCMGCCGCGGTAA; R=GGACTACHVGGGTWTCTAAT ;;
    V3V4) NAME=341F-805R; F=CCTACGGGNGGCWGCAG;   R=GACTACHVGGGTATCTAATCC ;;
    *) echo "REGION must be V4 or V3V4" >&2; exit 1 ;;
esac
log "extract-reads $NAME"
qiime feature-classifier extract-reads \
    --i-sequences "$ASSETS/gtdb-r214-ssu-seqs.qza" \
    --p-f-primer "$F" --p-r-primer "$R" \
    --p-n-jobs "$THREADS" \
    --o-reads "$ASSETS/gtdb-r214-$NAME-seqs.qza"
log "dereplicate (uniq)"
qiime rescript dereplicate \
    --i-sequences "$ASSETS/gtdb-r214-$NAME-seqs.qza" \
    --i-taxa "$ASSETS/gtdb-r214-ssu-tax.qza" --p-mode uniq --p-threads "$THREADS" \
    --o-dereplicated-sequences "$ASSETS/gtdb-r214-$NAME-seqs-derep.qza" \
    --o-dereplicated-taxa "$ASSETS/gtdb-r214-$NAME-tax-derep.qza"
log "fit-classifier-naive-bayes"
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$ASSETS/gtdb-r214-$NAME-seqs-derep.qza" \
    --i-reference-taxonomy "$ASSETS/gtdb-r214-$NAME-tax-derep.qza" \
    --o-classifier "$ASSETS/gtdb-r214-$NAME-classifier.qza"
log "done"
