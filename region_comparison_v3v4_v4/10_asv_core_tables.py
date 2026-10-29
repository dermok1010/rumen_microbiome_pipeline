"""ASV-level V4-core count tables for the dairy prediction test (11_dairy_kernel_asv.py).

Uses the 16S vs RE-RRS paper's cross-primer rule from micro_prediction_NZ
(re_rrs_comparison/build_nz_inputs.py, loco_16s_v4core): the V4 core is the first 250 bases after
the 515F primer, identical cores are summed and named by MD5.
  V3-V4 ASVs: cut after the 515F primer found by sequence (the paper's nz_core, imported).
  V4 ASVs:    the paper's ie_core assumes V4 ASVs start right after 515F, which holds for the Irish
              sheep but not for every dairy flow cell: CRT23 and part of the HM24/clover22 ASVs keep
              8-12 bases of the primer's 3' end. Those bases are removed first (the primer tail
              GCGGT[AC]A searched for in the first 20 bases), then the first 250 bases are kept.
Writes, to outputs/region_comparison_v3v4_v4/asv_core/:
  train_v4.tsv   HM24 + clover22 V4 samples (07's training samples are a subset)
  crt23_v4.tsv   CRT23 T1 and T2 V4 samples
  bsi_v3v4.tsv   BSI's Moorepark V3-V4 libraries
  core_stats.csv ASVs and reads kept per table
Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python (needs h5py).
"""
import hashlib
import re
import sys
import zipfile
from importlib import import_module
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(Path.home() / "microbiome_prediction" / "re_rrs_comparison")]
gt = import_module("02_genus_tables")
nzi = import_module("build_nz_inputs")
OUT = gt.OUT / "asv_core"
PRIMER_TAIL = re.compile("GCGGT[AC]A")


def v4_core(s):
    m = PRIMER_TAIL.search(s[:20])
    s = s[m.end():] if m else s
    return s[:nzi.CORE_LEN] if len(s) >= nzi.CORE_LEN else None


def fasta(qza):
    with zipfile.ZipFile(qza) as z:
        text = z.read(next(n for n in z.namelist() if n.endswith("dna-sequences.fasta"))).decode()
    return {r.split("\n", 1)[0]: r.split("\n", 1)[1].replace("\n", "") for r in text.split(">")[1:]}


def core_table(t, seqs, core):
    """t: ASVs x samples. Returns samples x cores (MD5 ids) and kept-ASV / kept-read shares."""
    ids = {a: core(seqs[a]) for a in t.index}
    keep = [a for a in t.index if ids[a]]
    out = t.loc[keep].groupby(lambda a: hashlib.md5(ids[a].encode()).hexdigest()).sum().T
    return out, len(keep) / len(t), t.loc[keep].to_numpy().sum() / t.to_numpy().sum()


def main():
    OUT.mkdir(exist_ok=True)
    cols = pd.read_csv(gt.V4_TABLE, sep="\t", skiprows=1, nrows=0).columns
    sel = [c for c in cols if c.startswith(("EN00011687__Dairy_Trt", "clover22__", "CRT23__"))]
    v4 = pd.read_csv(gt.V4_TABLE, sep="\t", skiprows=1, usecols=[cols[0]] + sel, index_col=0)
    v4.index.name = "feature-id"
    v4 = v4[(v4 > 0).any(axis=1)]
    v4_seqs = fasta(nzi.IE_SEQS)
    v34 = gt.v3v4_table()
    v34_seqs = fasta(nzi.NZ_SEQS)
    stats = []
    for name, t, seqs, core in [
        ("train_v4", v4[[c for c in v4.columns if not c.startswith("CRT23__")]], v4_seqs, v4_core),
        ("crt23_v4", v4[[c for c in v4.columns if c.startswith("CRT23__")]], v4_seqs, v4_core),
        ("bsi_v3v4", v34, v34_seqs, nzi.nz_core)]:
        t = t[(t > 0).any(axis=1)]
        tab, asvs_kept, reads_kept = core_table(t, seqs, core)
        tab.index.name = "sample"
        tab.reset_index().to_csv(OUT / f"{name}.tsv", sep="\t", index=False)
        stats.append({"table": name, "samples": len(tab), "asvs": len(t), "cores": tab.shape[1],
                      "pct_asvs_kept": round(100 * asvs_kept, 2), "pct_reads_kept": round(100 * reads_kept, 2)})
    stats = pd.DataFrame(stats)
    stats.to_csv(OUT / "core_stats.csv", index=False)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
