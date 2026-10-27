#!/bin/bash
# GTDB r214 (V4 classifier) taxonomy for the V4 ASVs of the other dairy trials (HM24, clover22),
# the training herds for the prediction test (07_dairy_prediction.py). ASVs already classified for
# CRT23 by 01_classify.sh are skipped. Run from the repo root.
source "$(dirname "$0")/../classifier/build/gtdb_r214/common.sh"
OUT=$REPO/outputs/region_comparison_v3v4_v4
"$PY" - "$OUT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
out = Path(sys.argv[1])
f = "results/rumen_combined_feature-table.tsv"
cols = pd.read_csv(f, sep="\t", skiprows=1, nrows=0).columns
sel = [c for c in cols if c.startswith(("EN00011687__Dairy_Trt", "clover22__"))]
t = pd.read_csv(f, sep="\t", skiprows=1, usecols=[cols[0]] + sel, index_col=0)
done = pd.read_csv(out / "taxonomy_V4_gtdb_r214.tsv", sep="\t", index_col=0).index
need = t.index[(t > 0).any(axis=1) & ~t.index.isin(done)]
pd.DataFrame(index=pd.Index(need, name="feature-id")).to_csv(out / "asvs_V4_dairy_train.tsv", sep="\t")
print(len(sel), "samples,", len(need), "ASVs to classify")
PY
qiime feature-table filter-seqs --i-data "$REPO/outputs/rumen_combined_rep_seqs.qza" \
    --m-metadata-file "$OUT/asvs_V4_dairy_train.tsv" --o-filtered-data "$OUT/rep_seqs_V4_dairy_train.qza"
log "classify-sklearn with gtdb-r214-515F-806R"
qiime feature-classifier classify-sklearn \
    --i-classifier "$ASSETS/gtdb-r214-515F-806R-classifier.qza" --i-reads "$OUT/rep_seqs_V4_dairy_train.qza" \
    --p-n-jobs "$THREADS" --p-reads-per-batch 2000 --p-read-orientation same \
    --o-classification "$OUT/taxonomy_V4_dairy_train_gtdb_r214.qza"
qiime tools export --input-path "$OUT/taxonomy_V4_dairy_train_gtdb_r214.qza" --output-path "$OUT/exp_train"
mv "$OUT/exp_train/taxonomy.tsv" "$OUT/taxonomy_V4_dairy_train_gtdb_r214.tsv" && rmdir "$OUT/exp_train"
log "done"
