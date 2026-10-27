"""Genus tables (GTDB r214) for the V3-V4 vs V4 dairy comparison.

  --write-ids  write the ASV lists to classify (V4: CRT23 samples; V3-V4: the Moorepark dairy
               samples in BSI's comparison table), for 01_classify.sh
  (default)    collapse both regions' ASV tables to GTDB genus, one row per sample, plus a
               sample sheet linking every sample to its cow

Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python.
"""
import argparse
import io
import re
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "region_comparison_v3v4_v4"
V4_TABLE = REPO / "results" / "rumen_combined_feature-table.tsv"
CMP = Path("/home/dermodkkelly/microbiome_prediction/data/re_rrs/comparison_16s_vs_rerrs_2024")
V3V4_TABLE = CMP / "16S_table" / "table.qza"  # BSI samples, Dermot's HPC run 2026-08-30
UNASSIGNED = "Unassigned_genus"


def v4_table():
    cols = pd.read_csv(V4_TABLE, sep="\t", skiprows=1, nrows=0).columns
    crt = [c for c in cols if c.startswith("CRT23__")]
    t = pd.read_csv(V4_TABLE, sep="\t", skiprows=1, usecols=[cols[0]] + crt, index_col=0)
    t.index.name = "feature-id"
    return t[(t > 0).any(axis=1)]


def dairy_samples():
    """BSI comparison samples of Moorepark dairy cows: 16S sample id -> cow number."""
    info = pd.read_csv(CMP / "sampleInfo.csv")
    dairy = info[(info.Species == "CATTLE") & info.COMMENT.astype(str).str.startswith("Moorepark")]
    cow = dict(zip(dairy.LabSampleID.astype(str), dairy.UIDTag.astype(str)))
    return cow


def v3v4_table():
    z = zipfile.ZipFile(V3V4_TABLE)
    name = next(n for n in z.namelist() if n.endswith("feature-table.biom"))
    with h5py.File(io.BytesIO(z.read(name)), "r") as f:
        obs = [x.decode() for x in f["observation/ids"][:]]
        smp = [x.decode() for x in f["sample/ids"][:]]
        data, ind, ptr = f["sample/matrix/data"][:], f["sample/matrix/indices"][:], f["sample/matrix/indptr"][:]
    cow = dairy_samples()
    # replicate library preps carry a trailing letter (946397a); keep them, mapped to the same cow
    keep = [i for i, s in enumerate(smp) if re.sub(r"[a-z]+$", "", s) in cow]
    rows = {}
    for i in keep:  # CSR by sample: row i's counts are data[ptr[i]:ptr[i+1]] at columns ind[...]
        rows[smp[i]] = pd.Series(data[ptr[i]:ptr[i + 1]], index=[obs[j] for j in ind[ptr[i]:ptr[i + 1]]])
    t = pd.DataFrame(rows).fillna(0)
    t.index.name = "feature-id"
    return t


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
    tables = {"V4": v4_table(), "V3V4": v3v4_table()}

    if args.write_ids:
        for region, t in tables.items():
            pd.DataFrame(index=t.index).to_csv(OUT / f"asvs_{region}.tsv", sep="\t")
            print(f"{region}: {t.shape[1]} samples, {len(t)} ASVs")
        return

    cow = dairy_samples()
    sheet = []
    for region, t in tables.items():
        tax = pd.read_csv(OUT / f"taxonomy_{region}_gtdb_r214.tsv", sep="\t", index_col=0)["Taxon"]
        genus = t.index.map(lambda a: genus_of(tax.get(a, "")))
        g = t.groupby(genus).sum().T
        g.index.name = "sample"
        g.reset_index().to_csv(OUT / f"genus_{region}_gtdb_r214.tsv", sep="\t", index=False)
        for s in g.index:
            if region == "V4":
                tp, c = re.match(r"CRT23__(T\d)\.(\d+)", s).groups()
            else:
                tp, c = "BSI", cow[re.sub(r"[a-z]+$", "", s)]
            sheet.append({"sample": s, "region": region, "timepoint": tp, "cow": c,
                          "replicate": region == "V3V4" and s[-1].isalpha(),
                          "reads": int(g.loc[s].sum()),
                          "pct_to_genus": round(100 * (1 - g.loc[s].get(UNASSIGNED, 0) / g.loc[s].sum()), 2)})
        print(region, g.shape)
    pd.DataFrame(sheet).to_csv(OUT / "samples.csv", index=False)


if __name__ == "__main__":
    main()
