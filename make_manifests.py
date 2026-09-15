#!/usr/bin/env python3
"""
Generate per-flowcell QIIME2 manifests (PairedEndFastqManifestPhred33V2 / tab-separated)
for a combined multi-batch run, or a deliberate subset of it.

- Groups files by SEQUENCING RUN (instrument:run:flowcell read from the FASTQ header),
  so DADA2 can denoise each run separately (correct error model per run).
- sample-id = <BATCH>__<filename-stem>, guaranteeing global uniqueness and
  encoding the batch as a covariate. The original stem is preserved for
  later metadata joins (just strip the '<BATCH>__' prefix).
- Pairs forward (_1) and reverse (_2) by the stem before the FINAL _1/_2.
- ERRORS OUT if any forward lacks a reverse partner (or vice versa).

Outputs:
  <outdir>/manifest_<RUNKEY>.tsv        one per flowcell
  <outdir>/sampleid_to_file_map.tsv     full audit: sample-id -> batch, stem, R1, R2, run
  <outdir>/run_summary.tsv              per-run file counts

Subsetting (for VM/fixture use -- "same/different/some of the same data"
without a manual manifest edit):
  --batch EN00011687 [--batch ...]      only these batch labels
  --flowcell M07726_229 [...]           only these flowcell run keys (post
                                         hoc filter, applied after reading
                                         headers -- still has to touch every
                                         file's header once to know)
  --sample-id-file path.txt             only sample-ids listed in this file
                                         (one per line, matching the
                                         <BATCH>__<stem> convention)
"""
import argparse
import gzip
import glob
import os
import re
import sys
from collections import defaultdict

DEFAULT_DATA_ROOT = os.environ.get(
    "RUMEN_PIPELINE_RAW_DATA_ROOT", "/data/Genetics/primary/R1240_microbiome"
)

NZ_METHANE_ROOT = "/data/BioScience/primary/R2002_methanepredict/NZ_method_comparison"

# batch folder -> batch label used in the sample-id prefix, or
# (label, root) when the folder lives under a different root than
# --data-root (e.g. Tully_16S, which is on a different HPC data share
# entirely). This is the one piece of real domain knowledge that has to be
# hardcoded somewhere (someone has to know which folder is which cohort) --
# kept as an explicit, documented registry rather than inferred from folder
# names.
BATCHES = {
    "dairy_20260608":                  "EN00010710",
    "sheep_dairy_ct_2024_EN00011679":  "EN00011679",
    "sheep_ct_2024_EN00011681":        "EN00011681",
    "sheep_ct_2024_hoggets_EN00011682":"EN00011682",
    "sheep_inzac_2025_EN00011685":     "EN00011685",
    "beef_sheep_2023_2024_EN00011686": "EN00011686",
    "dairy_beef_2023_EN00011687":      "EN00011687",
    "sheep_ct_2024_etc_EN00011684":    "EN00011684",
    "dairy_20260629_EN00012132":       "EN00012132",
    "sheep_inzac_EN00011689":          "EN00011689",
    "CT_microbiome_data":              "CTmicro",
    # NZ_comparison (inside R1240_microbiome) deliberately excluded (different platform)
    "clover22_16SV4":                  "clover22",
    "CRT23_16SV4":                     "CRT23",
    "Tully_16S":                       ("Tully", NZ_METHANE_ROOT),
}


def flowcell_of(path):
    """Read first header line, return instrument:run:flowcell."""
    with gzip.open(path, "rt") as fh:
        h = fh.readline().strip()
    # @INSTR:RUN:FLOWCELL:LANE:...   -> take fields 1-3 (strip leading @)
    parts = h.split(":")
    if len(parts) < 3:
        return "UNKNOWN"
    instr = parts[0].lstrip("@")
    return f"{instr}_{parts[1]}_{parts[2]}"


# Macrogen delivery folders (clover22_16SV4, CRT23_16SV4) carry pre-trim
# duplicates (<stem>.raw_1/.raw_2.fastq.gz in 01.RawData/) and
# merged/overlapped single reads (<stem>.extendedFrags.fastq.gz in
# 01.RawData/, <stem>.effective.fastq.gz in 00.CleanData/) alongside the
# delivered trimmed pair (<stem>_1/_2.fastq.gz) in the same directory tree.
# Only the trimmed pair should feed the pipeline -- confirmed by diffing
# read counts (identical) and lengths (.raw is uniformly full-cycle-length,
# untrimmed) for one sample. Skip the others outright rather than letting
# them register as spurious extra samples or unpaired-file errors.
EXCLUDE_PATTERNS = (".raw_1.fastq.gz", ".raw_2.fastq.gz", ".extendedFrags.fastq.gz",
                    ".effective.fastq.gz")


def stem_and_read(fname):
    """From e.g. Sheep_CT24_39_1.fastq.gz -> ('Sheep_CT24_39', '1').
       Handles replicate names like 3097_1_1 -> stem '3097_1', read '1'.
       Also handles bcl2fastq naming, e.g.
       10251_S94_L001_R2_001.fastq.gz -> ('10251', '2')."""
    base = fname[:-len(".fastq.gz")] if fname.endswith(".fastq.gz") else fname
    m = re.match(r"^(.+)_S\d+_L\d+_R([12])_\d+$", base)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^(.*)_([12])$", base)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", nargs="?", default="./manifest_files_combined")
    ap.add_argument("--data-root", default=DEFAULT_DATA_ROOT,
                     help="Root directory containing the per-batch FASTQ folders "
                          "(default: $RUMEN_PIPELINE_RAW_DATA_ROOT or the HPC path)")
    ap.add_argument("--batch", action="append", default=None,
                     help="Only include this batch label (e.g. EN00011687). Repeatable.")
    ap.add_argument("--flowcell", action="append", default=None,
                     help="Only include this flowcell run key (e.g. M07726_229). Repeatable.")
    ap.add_argument("--sample-id-file", default=None,
                     help="Only include sample-ids listed in this file (one per line).")
    ap.add_argument("--known-orphan", action="append", default=[],
                     help="Sample-id (e.g. Tully__50817) confirmed to genuinely be missing "
                          "its mate in the source delivery, not a naming/pairing bug. Skips "
                          "it with a note instead of failing the whole build. Repeatable. "
                          "Don't use this to paper over an unexpected pairing failure -- "
                          "confirm the mate is actually absent from the source first.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    known_orphans = set(args.known_orphan)

    sample_id_filter = None
    if args.sample_id_file:
        with open(args.sample_id_file) as f:
            sample_id_filter = {line.strip() for line in f if line.strip()}

    def label_of(spec):
        return spec[0] if isinstance(spec, tuple) else spec

    batches = BATCHES
    if args.batch:
        wanted = set(args.batch)
        batches = {k: v for k, v in BATCHES.items() if label_of(v) in wanted}
        if not batches:
            sys.stderr.write(f"--batch matched nothing in BATCHES: {sorted(wanted)}\n")
            sys.exit(1)

    # gather: run -> list of (sample_id, R1path, R2path, batch, stem)
    runs = defaultdict(list)
    audit = []
    errors = []

    for folder, spec in batches.items():
        label, root = spec if isinstance(spec, tuple) else (spec, args.data_root)
        bdir = os.path.join(root, folder)
        if not os.path.isdir(bdir):
            errors.append(f"MISSING FOLDER: {bdir}")
            continue
        files = glob.glob(os.path.join(bdir, "**", "*.fastq.gz"), recursive=True)
        files = [f for f in files if not any(f.endswith(p) for p in EXCLUDE_PATTERNS)]
        byfwd = {}
        byrev = {}
        for f in files:
            stem, read = stem_and_read(os.path.basename(f))
            if stem is None:
                errors.append(f"UNPARSEABLE NAME: {f}")
                continue
            (byfwd if read == "1" else byrev)[stem] = f
        all_stems = set(byfwd) | set(byrev)
        for stem in sorted(all_stems):
            sample_id = f"{label}__{stem}"
            if sample_id_filter is not None and sample_id not in sample_id_filter:
                continue
            r1 = byfwd.get(stem); r2 = byrev.get(stem)
            if r1 is None or r2 is None:
                if sample_id in known_orphans:
                    print(f"SKIPPING known orphan {sample_id}: R1={r1} R2={r2}", file=sys.stderr)
                else:
                    errors.append(f"UNPAIRED in {label}: stem '{stem}' R1={r1} R2={r2}")
                continue
            run = flowcell_of(r1)
            if args.flowcell and run not in set(args.flowcell):
                continue
            runs[run].append((sample_id, r1, r2))
            audit.append((sample_id, label, stem, run, r1, r2))

    # fail loudly before writing anything
    if errors:
        sys.stderr.write("\n!!! PROBLEMS FOUND -- no manifests written:\n")
        for e in errors:
            sys.stderr.write("  " + e + "\n")
        sys.exit(1)

    if not runs:
        sys.stderr.write("No samples matched the requested filters -- no manifests written.\n")
        sys.exit(1)

    # write per-run manifests
    for run, rows in sorted(runs.items()):
        path = os.path.join(args.outdir, f"manifest_{run}.tsv")
        with open(path, "w") as out:
            out.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for sid, r1, r2 in sorted(rows):
                out.write(f"{sid}\t{r1}\t{r2}\n")

    # audit + summary
    with open(os.path.join(args.outdir, "sampleid_to_file_map.tsv"), "w") as out:
        out.write("sample-id\tbatch\toriginal-stem\trun\tR1\tR2\n")
        for row in sorted(audit):
            out.write("\t".join(row) + "\n")

    with open(os.path.join(args.outdir, "run_summary.tsv"), "w") as out:
        out.write("run\tn_samples\tmanifest\n")
        for run, rows in sorted(runs.items()):
            out.write(f"{run}\t{len(rows)}\tmanifest_{run}.tsv\n")

    print(f"Wrote {len(runs)} per-run manifests to {args.outdir}")
    print(f"{'RUN (instr_run_flowcell)':32} {'n_samples':>9}")
    print("-" * 44)
    total = 0
    for run, rows in sorted(runs.items()):
        print(f"{run:32} {len(rows):>9}")
        total += len(rows)
    print("-" * 44)
    print(f"{'TOTAL':32} {total:>9}")
    print(f"\nAudit map: {args.outdir}/sampleid_to_file_map.tsv")
    print(f"Run summary: {args.outdir}/run_summary.tsv")


if __name__ == "__main__":
    main()
