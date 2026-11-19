#!/bin/bash
# Download the Greengenes2 2024.09 backbone (V4 515F/806R extract, full length, taxonomy) from the
# official release FTP into $REF (default ~/refs/greengenes2_2024.09) and record md5s.
# Vourlaki et al. 2026 (ISME J, doi 10.1093/ismejo/wrag185) used "a customized GreenGenes2
# (release 2024.09)" trained for each dataset's primers; their Irish data are V4.
source "$(dirname "$0")/common.sh"
BASE=http://ftp.microbio.me/greengenes_release/2024.09
mkdir -p "$REF"; cd "$REF"
for f in 00README 00CHANGELOG 00LICENSE 2024.09.backbone.v4.fna.qza 2024.09.backbone.tax.qza \
         2024.09.backbone.full-length.fna.qza; do
    [ -s "$f" ] || { log "get $f"; curl -sS -fL -o "$f" "$BASE/$f"; }
done
md5sum 00* 2024.09.* > MD5SUMS
{ echo "Source: $BASE/ (downloaded $(date -u '+%F %T') UTC)"; cat MD5SUMS; } > SOURCE.txt
cat SOURCE.txt
