"""GTDB r214 genus table for the 92 NZ Flock 2638 lambs (V3-V4, BSI comparison run).

  --write-ids  write the lambs' ASV list for 13_nz_sheep_classify.sh
  (default)    collapse their ASV counts to GTDB genus: genus_NZ_sheep_gtdb_r214.tsv, one row per
               lamb (SampleID as in micro_prediction_NZ's prep_inputs/nz_16s_asv), same layout as
               outputs/gtdb_r214_reclassification/<species>_genus_counts_gtdb_r214.tsv

Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

genus_of = import_module("02_genus_tables").genus_of
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "region_comparison_v3v4_v4"
NZ = Path("/home/dermodkkelly/microbiome_prediction/re_rrs_comparison/prep_inputs/nz_16s_asv/counts.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-ids", action="store_true")
    args = ap.parse_args()
    t = pd.read_csv(NZ, sep="\t", index_col="SampleID")
    t = t.loc[:, t.sum() > 0]
    if args.write_ids:
        pd.DataFrame(index=pd.Index(t.columns, name="feature-id")).to_csv(OUT / "asvs_NZ_sheep.tsv", sep="\t")
        print(f"NZ sheep: {len(t)} lambs, {t.shape[1]} ASVs")
        return
    tax = pd.read_csv(OUT / "taxonomy_NZ_sheep_gtdb_r214.tsv", sep="\t", index_col=0)
    g = t.T.groupby(tax.Taxon.reindex(t.columns).map(genus_of)).sum().T
    g.index.name = "SampleID"
    g.to_csv(OUT / "genus_NZ_sheep_gtdb_r214.tsv", sep="\t")
    unas = g.get("Unassigned_genus", 0).sum() / g.values.sum()
    print(f"NZ sheep: {g.shape[1]} genera, {100 * (1 - unas):.1f}% of reads to genus")


if __name__ == "__main__":
    main()
