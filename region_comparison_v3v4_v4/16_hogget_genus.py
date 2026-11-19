"""GTDB r214 genus table for the 91 phenotyped INZAC Irish hoggets (V3-V4, BSI comparison run).
A copy of 14_nz_sheep_genus.py for the hoggets.

  --write-ids  write the hoggets' ASV list for 15_hogget_classify.sh
  (default)    collapse their ASV counts to GTDB genus: genus_INZAC_hoggets_gtdb_r214.tsv, one row per
               hogget (SampleID as in micro_prediction_NZ's prep_inputs/hogcv_16s_asv), same layout as
               genus_NZ_sheep_gtdb_r214.tsv

Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python. PIPELINE_DATA_REPO (default:
this repo) is the checkout whose outputs/ are used.
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

genus_of = import_module("02_genus_tables").genus_of
REPO = Path(os.environ.get("PIPELINE_DATA_REPO", Path(__file__).resolve().parents[1]))
OUT = REPO / "outputs" / "region_comparison_v3v4_v4"
HOG = Path("/home/dermodkkelly/microbiome_prediction/re_rrs_comparison/prep_inputs/hogcv_16s_asv/counts.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-ids", action="store_true")
    args = ap.parse_args()
    t = pd.read_csv(HOG, sep="\t", index_col="SampleID")
    t = t.loc[:, t.sum() > 0]
    if args.write_ids:
        pd.DataFrame(index=pd.Index(t.columns, name="feature-id")).to_csv(OUT / "asvs_INZAC_hoggets.tsv", sep="\t")
        print(f"hoggets: {len(t)} animals, {t.shape[1]} ASVs")
        return
    tax = pd.read_csv(OUT / "taxonomy_INZAC_hoggets_gtdb_r214.tsv", sep="\t", index_col=0)
    assert t.columns.isin(tax.index).all()
    g = t.T.groupby(tax.Taxon.reindex(t.columns).map(genus_of)).sum().T
    g.index.name = "SampleID"
    g.to_csv(OUT / "genus_INZAC_hoggets_gtdb_r214.tsv", sep="\t")
    unas = g.get("Unassigned_genus", 0).sum() / g.values.sum()
    print(f"hoggets: {g.shape[1]} genera, {100 * (1 - unas):.1f}% of reads to genus")


if __name__ == "__main__":
    main()
