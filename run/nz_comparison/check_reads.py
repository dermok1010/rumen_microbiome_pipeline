#!/usr/bin/env python3
"""
Pre-flight check for the NZ_comparison 16S run, before any Slurm job.

Answers the three things the DADA2/classifier settings depend on, from a
sample of reads (standard library only, safe on a login node, ~1 min):

  1. Which primer pair is on the reads: 515F/806R (V4, this pipeline's
     default) or 341F/805R (V3-V4, which the 2026-09 one-off QIIME2 run's
     trim_left 17/21 suggests).
  2. Read length (2x250 vs 2x300), which bounds trunc_len_f/trunc_len_r.
  3. Where median quality drops, to choose trunc_len_f/trunc_len_r.

Also reports file counts, R1/R2 pairing by name, and flowcell(s).

Usage: python3 check_reads.py [FOLDER] [--files N] [--reads N]
"""
import argparse
import collections
import glob
import gzip
import os
import re
import statistics

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
PRIMERS = {
    "515F (V4 fwd)":    ("R1", "GTGYCAGCMGCCGCGGTAA"),
    "341F (V3-V4 fwd)": ("R1", "CCTACGGGNGGCWGCAG"),
    "806R (V4 rev)":    ("R2", "GGACTACNVGGGTWTCTAAT"),
    "805R (V3-V4 rev)": ("R2", "GACTACHVGGGTATCTAATCC"),
}
QUAL_POSITIONS = (150, 200, 230, 250, 260, 270, 280, 290, 300)


def rx(primer):
    # Allow up to 8 leading bases (heterogeneity spacers / phasing).
    return re.compile(r"^[ACGTN]{0,8}" + "".join(IUPAC[b] for b in primer))


def read_type(path):
    name = os.path.basename(path)
    if re.search(r"(_R1_\d+|_1)\.fastq\.gz$", name):
        return "R1"
    if re.search(r"(_R2_\d+|_2)\.fastq\.gz$", name):
        return "R2"
    return "?"


def stem(path):
    return re.sub(r"(_R[12]_\d+|_[12])\.fastq\.gz$", "", os.path.basename(path))


def scan(path, n_reads):
    seqs, quals, header = [], [], None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            k = i % 4
            if k == 0 and header is None:
                header = line.strip()
            elif k == 1:
                seqs.append(line.strip())
            elif k == 3:
                quals.append(line.strip())
                if len(quals) >= n_reads:
                    break
    return header, seqs, quals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default="/data/Genetics/primary/R1240_microbiome/NZ_comparison")
    ap.add_argument("--files", type=int, default=12, help="R1/R2 pairs to sample")
    ap.add_argument("--reads", type=int, default=2000, help="reads per file")
    a = ap.parse_args()

    files = sorted(f for f in glob.glob(os.path.join(a.folder, "**", "*.fastq.gz"), recursive=True)
                   if not re.search(r"\.(extendedFrags|effective)\.fastq\.gz$", f)
                   and "Undetermined" not in f)
    print(f"Folder: {a.folder}\nFASTQ files: {len(files)}")
    if not files:
        print("No FASTQs found -- check the folder path before going further.")
        return
    by = collections.defaultdict(dict)
    for f in files:
        by[stem(f)][read_type(f)] = f
    pairs = sorted(s for s, d in by.items() if "R1" in d and "R2" in d)
    unpaired = sorted(s for s, d in by.items() if not ("R1" in d and "R2" in d))
    print(f"Samples with R1+R2: {len(pairs)}   unpaired/unrecognised: {len(unpaired)}")
    for s in unpaired[:10]:
        print(f"  unpaired: {s} -> {sorted(by[s])}")
    print("Example names:", ", ".join(os.path.basename(by[s]['R1']) for s in pairs[:3]))

    step = max(1, len(pairs) // a.files)
    chosen = pairs[::step][:a.files]
    hits = collections.Counter()
    total = collections.Counter()
    lengths = {"R1": collections.Counter(), "R2": collections.Counter()}
    qual_at = {"R1": collections.defaultdict(list), "R2": collections.defaultdict(list)}
    flowcells = collections.Counter()
    for s in chosen:
        for rt in ("R1", "R2"):
            header, seqs, quals = scan(by[s][rt], a.reads)
            if header:
                p = header.lstrip("@").split(":")
                flowcells[":".join(p[:3]) if len(p) >= 3 else header] += 1
            total[rt] += len(seqs)
            for name, (prt, primer) in PRIMERS.items():
                if prt == rt:
                    r = rx(primer)
                    hits[name] += sum(1 for q in seqs if r.match(q))
            for q in seqs:
                lengths[rt][len(q)] += 1
            for pos in QUAL_POSITIONS:
                vals = [ord(q[pos - 1]) - 33 for q in quals if len(q) >= pos]
                if vals:
                    qual_at[rt][pos].append(statistics.median(vals))

    print(f"\nSampled {len(chosen)} pairs x {a.reads} reads")
    print("Flowcells (instrument:run:flowcell):", dict(flowcells))
    print("\nPrimer at read start (fraction of reads):")
    for name, (rt, _) in PRIMERS.items():
        frac = hits[name] / total[rt] if total[rt] else 0
        print(f"  {name:18s} {rt}: {frac:6.1%}")
    for rt in ("R1", "R2"):
        common = lengths[rt].most_common(3)
        print(f"\n{rt} read lengths (top 3): {common}")
        print(f"{rt} median Phred by position (median across sampled files):")
        print("   " + "  ".join(f"{pos}:{statistics.median(v):.0f}" for pos, v in sorted(qual_at[rt].items())))
    print("\nNext: fill run/nz_comparison/nz_comparison.params.yaml from this "
          "(see run/nz_comparison/README.md, step 2).")


if __name__ == "__main__":
    main()
