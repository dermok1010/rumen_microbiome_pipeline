"""Genus count tables from the GTDB r214 classification of the comparison ASVs.

  --write-ids  write needed_asvs.tsv: ASVs with any reads in the comparison animals (step 3 input)
  (default)    collapse the ASV table to GTDB genus for sheep and beef, one table each

The comparison animals are the 16S-vs-RE-RRS prediction population in micro_prediction_NZ
(prep_inputs/<species>_16s_genus/metadata.csv with a methane record: 505 sheep, 116 beef).
Genus labels are GTDB names without the "g__" prefix, matching BSI's RE-RRS genus tables.
ASVs not classified to genus (or missing from the rep seqs) are pooled as "Unassigned_genus".

Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python.
"""
import argparse
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "outputs" / "gtdb_r214_reclassification"
TABLE = REPO / "results" / "rumen_combined_feature-table.tsv"  # 1 Oct 2026 run
PREP = Path("/home/dermodkkelly/microbiome_prediction/re_rrs_comparison/prep_inputs")
TARGET = {"sheep": "ch4_g_day2_1v3", "beef": "CH4__g_d_"}
UNASSIGNED = "Unassigned_genus"


def comparison_samples():
    out = {}
    for species, target in TARGET.items():
        m = pd.read_csv(PREP / f"{species}_16s_genus" / "metadata.csv")
        out[species] = m.loc[m[target].notna(), "SampleID"].tolist()
    return out


def asv_table(sample_ids):
    cols = pd.read_csv(TABLE, sep="\t", skiprows=1, nrows=0).columns
    tab = pd.read_csv(TABLE, sep="\t", skiprows=1, usecols=[cols[0]] + sample_ids, index_col=0)
    tab.index.name = "feature-id"
    return tab[(tab > 0).any(axis=1)]


def genus_of(taxon):
    for rank in str(taxon).split(";"):
        rank = rank.strip()
        if rank.startswith("g__") and len(rank) > 3:
            return rank[3:]
    return UNASSIGNED


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-ids", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    samples = comparison_samples()
    tab = asv_table(samples["sheep"] + samples["beef"])

    if args.write_ids:
        pd.DataFrame(index=tab.index).to_csv(OUT / "needed_asvs.tsv", sep="\t")
        print(f"{len(tab)} ASVs in {tab.shape[1]} samples -> {OUT / 'needed_asvs.tsv'}")
        return

    tax = pd.read_csv(OUT / "comparison_taxonomy_gtdb_r214.tsv", sep="\t", index_col=0)
    genus = tab.index.map(lambda a: genus_of(tax["Taxon"].get(a, "")))
    summary = []
    for species, ids in samples.items():
        g = tab[ids].groupby(genus).sum().T
        g = g.loc[:, (g > 0).any(axis=0)]
        g.index.name = "SampleID"
        g.reset_index().to_csv(OUT / f"{species}_genus_counts_gtdb_r214.tsv", sep="\t", index=False)
        reads = g.to_numpy().sum()
        assigned = reads - g.get(UNASSIGNED, pd.Series(0)).sum()
        summary.append({"species": species, "samples": len(ids),
                        "genera": g.shape[1] - (UNASSIGNED in g.columns),
                        "pct_reads_to_genus": round(100 * assigned / reads, 2)})
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
