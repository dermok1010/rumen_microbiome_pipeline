"""Greengenes2 2024.09 genus tables for the combined 17-flowcell run (for the Ch5 Ruminosignatures
replication, ~/methane-signature-replication).

Same ASV table, control filtering and Unassigned pooling as gtdb_r214/06_full_genus_table.py and
postprocessing/01_build_microbiome_matrices.ipynb. Genus columns are the GG2 lineage down to genus
with rank prefixes stripped and GG2 numeric suffixes kept, joined by ';'
(e.g. "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales_737866;Succinivibrionaceae;UBA2810"),
which is exactly how Vourlaki et al. 2026 name the rows of their W matrices. ASVs with no genus
rank, an empty g__, or no taxonomy (no rep sequence) are pooled as "Unassigned".

Outputs (outputs/gg2_2024.09_full/):
  genus_counts_no_controls_gg2_2024.09.csv   SampleID x GG2 genus lineage (non-control samples)
  asv_taxonomy_gg2.tsv                       per ASV: reads, GG2 taxon, confidence, genus lineage
  summary.csv                                ASVs / reads assigned to genus
Run with ~/envs/rumen_microbiome_pipeline/bin/python from the repo root.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(os.environ.get("OUT", REPO / "outputs" / "gg2_2024.09_full"))
TABLE = REPO / "results" / "rumen_combined_feature-table.tsv"
UNASSIGNED = "Unassigned"
CONTROL = r'__(?:ControlN_ive|ControlNive|ControlP_ive|Minus|Plus|N|P)_?\d*$'
TULLY_CONTROL = r'^Tully__(?:16S|Zymo|N-|n-|NEG|H2O|saliva)'


def genus_lineage(taxon):
    if pd.isna(taxon):
        return UNASSIGNED
    ranks = [p.strip() for p in str(taxon).split(";")]
    if len(ranks) < 6 or not ranks[5].startswith("g__") or not ranks[5][3:].strip():
        return UNASSIGNED
    return ";".join(r[3:].strip() for r in ranks[:6])


def main():
    tab = pd.read_csv(TABLE, sep="\t", skiprows=1, index_col=0, dtype={"#OTU ID": str}).astype(np.int64)
    ids = pd.Series(tab.columns)
    ctrl = (ids.str.contains(CONTROL, regex=True) | ids.str.contains(TULLY_CONTROL, regex=True)).to_numpy()
    bio = tab.loc[:, ~ctrl]
    bio = bio[(bio > 0).any(axis=1)]
    tx = pd.read_csv(OUT / "all_taxonomy_gg2_2024.09.tsv", sep="\t", index_col=0)
    g = pd.Series([genus_lineage(tx["Taxon"].get(a, np.nan)) for a in bio.index], index=bio.index)
    wide = bio.groupby(g.to_numpy()).sum().T
    wide.index.name = "SampleID"
    wide.reset_index().to_csv(OUT / "genus_counts_no_controls_gg2_2024.09.csv", index=False)
    reads = bio.sum(axis=1)
    per = pd.DataFrame({"reads": reads, "taxon": tx["Taxon"].reindex(bio.index),
                        "confidence": tx["Confidence"].reindex(bio.index), "genus": g})
    per.to_csv(OUT / "asv_taxonomy_gg2.tsv", sep="\t")
    asg = g != UNASSIGNED
    s = pd.DataFrame([dict(asvs=len(g), asvs_with_seq=int(per.taxon.notna().sum()),
                           asvs_to_genus=int(asg.sum()), share_asvs_to_genus=asg.mean(),
                           share_reads_to_genus=reads[asg].sum() / reads.sum(),
                           n_genera=int(wide.shape[1] - (UNASSIGNED in wide.columns)),
                           samples=wide.shape[0])])
    s.to_csv(OUT / "summary.csv", index=False)
    print(s.T.to_string())
    for k in ["UBA2810", "Kandleria"]:
        hit = [c for c in wide.columns if c.endswith(";" + k)]
        print(k, hit, int(wide[hit].to_numpy().sum()) if hit else 0, "reads")


if __name__ == "__main__":
    main()
