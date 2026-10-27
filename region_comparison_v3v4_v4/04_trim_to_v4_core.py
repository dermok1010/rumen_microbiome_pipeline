"""Trim both regions' ASVs to the same V4 core, so they can be compared sequence for sequence.

The V3-V4 amplicon (341F-805R) contains the V4 region. The core kept here runs from just after the
515F primer to just before the 805R/806R site:
  V3-V4 ASVs: cut after the 515F primer (found by sequence; IUPAC-aware); they already end where
              805R begins.
  V4 ASVs:    CRT23 reads keep 8 bases of 515F at the start (CGCGGTAA) and ~10 of 806R at the
              end (GGATTAGA.A); both are removed.
Identical cores are merged (counts summed) and named by MD5 of the sequence, QIIME2's own rule, so
a core seen in both regions gets one ID. Classification of the cores (05_core.sh) uses the GTDB
r214 V4 classifier for both, so identical sequences get identical genera: any remaining difference
comes from the PCR primers or the lab/sequencing, not from which part of the gene was classified.

  (default)  write core_seqs.fasta, core_table_{V4,V3V4}.tsv and core_overlap.csv
  --genus    after 05_core.sh has classified the cores: genus_{V4,V3V4}_core_gtdb_r214.tsv
Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python.
"""
import argparse
import hashlib
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module  # noqa: E402

gt = import_module("02_genus_tables")
OUT = gt.OUT
IUPAC = {"N": "[ACGT]", "M": "[AC]", "W": "[AT]"}
F515 = re.compile("".join(IUPAC.get(c, c) for c in "GTGCCAGCMGCCGCGGTAA"))
V4_START = "CGCGGTAA"                      # last 8 bases of 515F left on CRT23 reads
V4_END = re.compile(r"GGATTAGA[ACGT]{1,4}$")  # start of the 806R site left on CRT23 reads


def fasta(qza):
    z = zipfile.ZipFile(qza)
    text = z.read(next(n for n in z.namelist() if n.endswith("dna-sequences.fasta"))).decode()
    return {r.split("\n", 1)[0]: r.split("\n", 1)[1].replace("\n", "") for r in text.split(">")[1:]}


def core_v3v4(s):
    m = F515.search(s)
    return s[m.end():] if m else None


def core_v4(s):
    if not s.startswith(V4_START):
        return None
    s = s[len(V4_START):]
    m = V4_END.search(s)
    return s[:m.start()] if m else None


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()


def collapse(table, cores):
    t = table.loc[table.index.isin(cores.keys())]
    return t.groupby(t.index.map(lambda a: md5(cores[a]))).sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", action="store_true")
    args = ap.parse_args()
    tables = {"V4": gt.v4_table(), "V3V4": gt.v3v4_table()}

    if args.genus:
        tax = pd.read_csv(OUT / "core_taxonomy_gtdb_r214.tsv", sep="\t", index_col=0)["Taxon"]
        for region in tables:
            t = pd.read_csv(OUT / f"core_table_{region}.tsv", sep="\t", index_col=0)
            g = t.groupby(t.index.map(lambda a: gt.genus_of(tax.get(a, "")))).sum().T
            g.index.name = "sample"
            g.reset_index().to_csv(OUT / f"genus_{region}_core_gtdb_r214.tsv", sep="\t", index=False)
            print(region, g.shape)
        return

    seqs = {"V4": fasta(OUT / "rep_seqs_V4.qza"), "V3V4": fasta(OUT / "rep_seqs_V3V4.qza")}
    fn = {"V4": core_v4, "V3V4": core_v3v4}
    stats, cores_all, collapsed = [], {}, {}
    for region in tables:
        t = tables[region]
        cores = {a: fn[region](s) for a, s in seqs[region].items()}
        cores = {a: c for a, c in cores.items() if c}
        reads_kept = t.loc[t.index.isin(cores.keys())].to_numpy().sum() / t.to_numpy().sum()
        collapsed[region] = collapse(t, cores)
        for c in cores.values():
            cores_all[md5(c)] = c
        lens = pd.Series([len(c) for c in cores.values()])
        stats.append({"region": region, "asvs": len(seqs[region]), "asvs_trimmed": len(cores),
                      "pct_reads_kept": round(100 * reads_kept, 2), "unique_cores": len(collapsed[region]),
                      "core_len_median": int(lens.median()), "core_len_p05": int(lens.quantile(.05)),
                      "core_len_p95": int(lens.quantile(.95))})
    a, b = collapsed["V4"], collapsed["V3V4"]
    both = a.index.intersection(b.index)
    for s in stats:
        t = a if s["region"] == "V4" else b
        s["cores_also_in_other_region"] = len(both)
        s["pct_reads_in_shared_cores"] = round(100 * t.loc[both].to_numpy().sum() / t.to_numpy().sum(), 2)
    pd.DataFrame(stats).to_csv(OUT / "core_overlap.csv", index=False)
    print(pd.DataFrame(stats).to_string(index=False))
    for region, t in collapsed.items():
        t.index.name = "feature-id"
        t.to_csv(OUT / f"core_table_{region}.tsv", sep="\t")
    with open(OUT / "core_seqs.fasta", "w") as f:
        for h, s in cores_all.items():
            f.write(f">{h}\n{s}\n")


if __name__ == "__main__":
    main()
