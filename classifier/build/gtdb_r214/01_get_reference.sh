#!/bin/bash
# Download GTDB 214.1 SSU (species representatives, Bacteria + Archaea) as QIIME2 artifacts.
source "$(dirname "$0")/common.sh"
log "get-gtdb-data 214.1 SpeciesReps Both"
qiime rescript get-gtdb-data --p-version 214.1 --p-domain Both --p-db-type SpeciesReps \
    --o-gtdb-taxonomy "$ASSETS/gtdb-r214-ssu-tax.qza" \
    --o-gtdb-sequences "$ASSETS/gtdb-r214-ssu-seqs.qza" --verbose
log "done"
