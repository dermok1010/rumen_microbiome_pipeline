"""Full GTDB r214 genus table for the combined 17-flowcell run, plus SILVA reconciliation and a
SILVA->GTDB crosswalk (for rumen-core's move to GTDB r214).

  --write-ids  needed_asvs.tsv (ASVs with reads in a non-control sample AND a rep sequence; input of
               05_classify_all_asvs.sh) and asvs_without_rep_seqs.tsv (with reads, no sequence)
  (default)    genus tables, reconciliation, crosswalk (needs all_taxonomy_gtdb_r214.tsv)

Source data (same as postprocessing/01_build_microbiome_matrices.ipynb, which produced the SILVA
genus table rumen-core uses, postprocessing/output/genus_counts_no_controls_515F_806R.csv):
  results/rumen_combined_feature-table.tsv      ASV table, 1 Oct 2026 run (1969 samples)
  results/rumen_combined_exported-taxonomy.tsv  SILVA 138 taxonomy of the same ASVs
  outputs/rumen_combined_rep_seqs.qza           rep seqs, 15 Sep run (the 1 Oct rep seqs stayed on
                                                the HPC; feature IDs are sequence hashes, so the Sep
                                                file covers every ASV the two runs share)

Genus handling mirrors notebook 01: the genus label is the g__ rank with the prefix removed; ASVs with
no genus rank (classified only to family or above), with an empty g__, or with no taxonomy at all
(here: no rep sequence) are pooled into one "Unassigned" column. GTDB suffixes (_A, _B, ...) are kept
as distinct genera. Control libraries are removed with notebook 01's two regular expressions.

Outputs (outputs/gtdb_r214_full/):
  genus_counts_with_controls_gtdb_r214.csv / genus_counts_no_controls_gtdb_r214.csv
  asv_taxonomy_silva_gtdb.tsv        per ASV: reads, SILVA and GTDB taxonomy and genus
  silva_gtdb_genus_crosswalk.csv     per SILVA genus x GTDB genus: reads, ASVs, shares (both ways)
  reconciliation.csv / reconciliation.md
Run with ~/envs/rumen_microbiome_pipeline/bin/python from the repo root.
"""
import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(os.environ.get("OUT", REPO / "outputs" / "gtdb_r214_full"))
TABLE = REPO / "results" / "rumen_combined_feature-table.tsv"
SILVA_TAX = REPO / "results" / "rumen_combined_exported-taxonomy.tsv"
SILVA_GENUS = REPO / "postprocessing" / "output" / "genus_counts_no_controls_515F_806R.csv"
RC_SUBSET = Path(os.environ.get("RUMEN_CORE_ROOT", Path.home() / "rumen-core")) / \
    "data" / "phenotyped_subset_v2" / "genus_counts_phenotyped_subset.csv"
UNASSIGNED = "Unassigned"
# notebook 01's control patterns, verbatim
CONTROL = r'__(?:ControlN_ive|ControlNive|ControlP_ive|Minus|Plus|N|P)_?\d*$'
TULLY_CONTROL = r'^Tully__(?:16S|Zymo|N-|n-|NEG|H2O|saliva)'


def read_table():
    tab = pd.read_csv(TABLE, sep="\t", skiprows=1, index_col=0, dtype={"#OTU ID": str})
    tab.index.name = "FeatureID"
    tab = tab.astype(np.int64)
    return tab


def is_control(ids):
    s = pd.Series(ids)
    return (s.str.contains(CONTROL, regex=True) | s.str.contains(TULLY_CONTROL, regex=True)).to_numpy()


def fasta_ids(path):
    with open(path) as f:
        return {line[1:].strip().split()[0] for line in f if line.startswith(">")}


def genus_of(taxon):
    """notebook 01's parse_taxonomy, Genus rank only."""
    if pd.isna(taxon):
        return UNASSIGNED
    g = UNASSIGNED
    for p in str(taxon).split(";"):
        p = p.strip()
        if p.startswith("g__"):
            v = p[3:].strip()
            g = v if v else UNASSIGNED
    return g


def rank_of(taxon, prefix):
    if pd.isna(taxon):
        return ""
    for p in str(taxon).split(";"):
        p = p.strip()
        if p.startswith(prefix):
            return p[len(prefix):].strip()
    return ""


def collapse(tab, genus):
    g = tab.groupby(genus.reindex(tab.index).fillna(UNASSIGNED).to_numpy()).sum()
    w = g.T
    w.index.name = "SampleID"
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-ids", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tab = read_table()
    ctrl = is_control(tab.columns)
    bio = tab.loc[:, ~ctrl]
    with_reads = bio.index[(bio > 0).any(axis=1)]
    have_seq = fasta_ids(OUT / "rep_seqs_sep15.fasta")
    print(f"ASV table: {tab.shape[0]} ASVs x {tab.shape[1]} samples ({ctrl.sum()} controls); "
          f"{len(with_reads)} ASVs with reads in non-control samples; {len(have_seq)} rep seqs")

    if args.write_ids:
        need = [a for a in with_reads if a in have_seq]
        miss = [a for a in with_reads if a not in have_seq]
        pd.DataFrame({"feature-id": need}).to_csv(OUT / "needed_asvs.tsv", sep="\t", index=False)
        pd.DataFrame({"feature-id": miss, "reads": bio.loc[miss].sum(axis=1).to_numpy()}).to_csv(
            OUT / "asvs_without_rep_seqs.tsv", sep="\t", index=False)
        print(f"{len(need)} ASVs to classify; {len(miss)} without a rep sequence")
        return

    # ---- taxonomies ----
    gt = pd.read_csv(OUT / "all_taxonomy_gtdb_r214.tsv", sep="\t", index_col=0)["Taxon"]
    sv = pd.read_csv(SILVA_TAX, sep="\t", index_col=0)["Taxon"]
    g_gtdb = pd.Series([genus_of(gt.get(a, np.nan)) for a in tab.index], index=tab.index)
    g_silva = pd.Series([genus_of(sv.get(a, np.nan)) for a in tab.index], index=tab.index)

    # ---- GTDB genus tables (same layout as notebook 01) ----
    wide = collapse(tab, g_gtdb)
    wide.reset_index().to_csv(OUT / "genus_counts_with_controls_gtdb_r214.csv", index=False)
    wide_bio = wide.loc[~ctrl]
    wide_bio.reset_index().to_csv(OUT / "genus_counts_no_controls_gtdb_r214.csv", index=False)

    # ---- reconciliation against the SILVA genus table rumen-core uses ----
    silva_tab = pd.read_csv(SILVA_GENUS, index_col=0)
    silva_rebuilt = collapse(bio, g_silva)
    rc_ids = pd.read_csv(RC_SUBSET, usecols=["SampleID", "species"])
    rows = []
    same_ids = set(silva_tab.index) == set(wide_bio.index)
    rebuilt = silva_rebuilt.reindex(index=silva_tab.index, columns=silva_tab.columns).fillna(0)
    max_cell_diff = float((rebuilt - silva_tab).abs().to_numpy().max())
    extra_cols = sorted(set(silva_rebuilt.columns) - set(silva_tab.columns))
    tot_silva = silva_tab.sum(axis=1)
    tot_gtdb = wide_bio.reindex(silva_tab.index).sum(axis=1)
    tot_asv = bio.sum(axis=0).reindex(silva_tab.index)
    assigned_silva = tot_silva - silva_tab.get(UNASSIGNED, 0)
    assigned_gtdb = tot_gtdb - wide_bio.reindex(silva_tab.index).get(UNASSIGNED, 0)
    asv_assigned_silva = bio.loc[g_silva[g_silva != UNASSIGNED].index].sum(axis=0).reindex(silva_tab.index)
    miss = pd.read_csv(OUT / "asvs_without_rep_seqs.tsv", sep="\t")
    miss_reads = bio.loc[miss["feature-id"]].sum(axis=0).reindex(silva_tab.index)
    per_sample = pd.DataFrame({
        "reads_asv_table": tot_asv, "reads_silva_genus_table": tot_silva,
        "reads_gtdb_genus_table": tot_gtdb, "silva_genus_assigned_from_asv": asv_assigned_silva,
        "silva_genus_assigned_from_table": assigned_silva, "gtdb_genus_assigned": assigned_gtdb,
        "reads_asvs_without_rep_seq": miss_reads})
    per_sample["in_rumen_core_subset"] = per_sample.index.isin(rc_ids["SampleID"])
    per_sample = per_sample.join(rc_ids.set_index("SampleID")["species"])
    per_sample.to_csv(OUT / "reconciliation_per_sample.csv")

    def summ(d, label):
        return {"set": label, "samples": len(d),
                "max_abs_total_diff_asv_vs_silva": float((d.reads_asv_table - d.reads_silva_genus_table).abs().max()),
                "max_abs_total_diff_gtdb_vs_silva": float((d.reads_gtdb_genus_table - d.reads_silva_genus_table).abs().max()),
                "max_abs_assigned_diff_silva": float((d.silva_genus_assigned_from_asv - d.silva_genus_assigned_from_table).abs().max()),
                "total_reads": int(d.reads_asv_table.sum()),
                "pct_reads_genus_silva": round(100 * d.silva_genus_assigned_from_table.sum() / d.reads_asv_table.sum(), 2),
                "pct_reads_genus_gtdb": round(100 * d.gtdb_genus_assigned.sum() / d.reads_asv_table.sum(), 2),
                "pct_reads_no_rep_seq": round(100 * d.reads_asvs_without_rep_seq.sum() / d.reads_asv_table.sum(), 4),
                "max_pct_reads_no_rep_seq_one_sample": round(100 * (d.reads_asvs_without_rep_seq / d.reads_asv_table).max(), 3)}
    rec = [summ(per_sample, "all non-control samples")]
    sub = per_sample[per_sample.in_rumen_core_subset]
    rec.append(summ(sub, "rumen-core phenotyped subset v2"))
    for sp, d in sub.groupby("species"):
        rec.append(summ(d, f"subset: {sp}"))
    rec = pd.DataFrame(rec)
    rec.to_csv(OUT / "reconciliation.csv", index=False)

    n_gen_gtdb = int((wide_bio.columns != UNASSIGNED).sum())
    n_gen_silva = int((silva_tab.columns != UNASSIGNED).sum())
    with open(OUT / "reconciliation.md", "w") as f:
        f.write("# GTDB r214 full reclassification: reconciliation with the SILVA genus table\n\n")
        f.write(f"- ASV table: {TABLE.relative_to(REPO)} ({tab.shape[0]} ASVs, {tab.shape[1]} samples, "
                f"{int(ctrl.sum())} controls removed -> {int((~ctrl).sum())}).\n")
        f.write(f"- Same sample IDs as {SILVA_GENUS.relative_to(REPO)}: {same_ids} "
                f"({len(silva_tab)} samples).\n")
        f.write(f"- Rebuilding the SILVA genus table from this ASV table and "
                f"{SILVA_TAX.relative_to(REPO)} with notebook 01's rules: max |cell difference| = "
                f"{max_cell_diff:g}; genus columns not in the SILVA table: {len(extra_cols)}.\n")
        f.write(f"- ASVs with reads in non-control samples: {len(with_reads)}; without a rep sequence "
                f"(left Unassigned): {len(miss)}.\n")
        f.write(f"- Genera (excluding Unassigned): SILVA {n_gen_silva}, GTDB {n_gen_gtdb}.\n\n")
        f.write("| " + " | ".join(rec.columns) + " |\n|" + "---|" * rec.shape[1] + "\n")
        for _, r in rec.iterrows():
            f.write("| " + " | ".join(str(v) for v in r.values) + " |\n")
        f.write("\n")

    # ---- ASV-level taxonomy table and SILVA->GTDB crosswalk ----
    rc_bio = [s for s in bio.columns if s in set(rc_ids["SampleID"])]
    reads_all = bio.sum(axis=1)
    reads_rc = bio[rc_bio].sum(axis=1)
    asv = pd.DataFrame({
        "reads_noncontrol": reads_all, "reads_rumen_core_subset": reads_rc,
        "has_rep_seq": [a in have_seq for a in tab.index],
        "silva_taxon": [sv.get(a, "") for a in tab.index], "silva_genus": g_silva,
        "gtdb_taxon": [gt.get(a, "") for a in tab.index], "gtdb_genus": g_gtdb,
        "gtdb_family": [rank_of(gt.get(a, np.nan), "f__") for a in tab.index]})
    asv = asv[asv.reads_noncontrol > 0]
    asv.index.name = "FeatureID"
    asv.to_csv(OUT / "asv_taxonomy_silva_gtdb.tsv", sep="\t")

    cw = asv.groupby(["silva_genus", "gtdb_genus"]).agg(
        reads_rumen_core_subset=("reads_rumen_core_subset", "sum"),
        reads_noncontrol=("reads_noncontrol", "sum"),
        n_asvs=("reads_noncontrol", "size"),
        gtdb_family=("gtdb_family", lambda x: x.mode().iat[0] if len(x) else "")).reset_index()
    for col, base in [("share_of_silva_genus", "silva_genus"), ("share_of_gtdb_genus", "gtdb_genus")]:
        cw[col] = cw.reads_rumen_core_subset / cw.groupby(base).reads_rumen_core_subset.transform("sum")
    cw = cw.sort_values(["silva_genus", "reads_rumen_core_subset"], ascending=[True, False])
    cw.to_csv(OUT / "silva_gtdb_genus_crosswalk.csv", index=False)
    top = cw.loc[cw.groupby("silva_genus").reads_rumen_core_subset.idxmax()]
    top.to_csv(OUT / "silva_gtdb_genus_crosswalk_top.csv", index=False)
    print(rec.to_string(index=False))
    print(f"max cell diff SILVA rebuild: {max_cell_diff}; GTDB genera {n_gen_gtdb}, SILVA {n_gen_silva}")
    for g in ["Methanobrevibacter", "Shuttleworthia", "Kandleria", "Prevotella",
              "Christensenellaceae_R-7_group", "Methanosphaera"]:
        print(cw[cw.silva_genus == g].head(4)[["silva_genus", "gtdb_genus", "share_of_silva_genus"]]
              .to_string(index=False, header=False))


if __name__ == "__main__":
    main()
