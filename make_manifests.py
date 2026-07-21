#!/usr/bin/env python3
"""
Generate per-flowcell QIIME2 manifests (PairedEndFastqManifestPhred33V2 / tab-separated)
for a combined multi-batch run.

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
"""
import os, sys, gzip, glob, re
from collections import defaultdict

DATA_ROOT = "/data/Genetics/primary/R1240_microbiome"

# batch folder -> batch label used in the sample-id prefix
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
    # NZ_comparison deliberately excluded (different platform)
}

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "./manifest_files_combined"
os.makedirs(OUTDIR, exist_ok=True)

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

def stem_and_read(fname):
    """From e.g. Sheep_CT24_39_1.fastq.gz -> ('Sheep_CT24_39', '1').
       Handles replicate names like 3097_1_1 -> stem '3097_1', read '1'."""
    base = fname[:-len(".fastq.gz")] if fname.endswith(".fastq.gz") else fname
    m = re.match(r"^(.*)_([12])$", base)
    if not m:
        return None, None
    return m.group(1), m.group(2)

# gather: run -> list of (sample_id, R1path, R2path, batch, stem)
runs = defaultdict(list)
audit = []
errors = []

for folder, label in BATCHES.items():
    bdir = os.path.join(DATA_ROOT, folder)
    if not os.path.isdir(bdir):
        errors.append(f"MISSING FOLDER: {bdir}")
        continue
    files = glob.glob(os.path.join(bdir, "**", "*.fastq.gz"), recursive=True)
    # index by stem
    byfwd = {}   # stem -> R1 path
    byrev = {}   # stem -> R2 path
    for f in files:
        stem, read = stem_and_read(os.path.basename(f))
        if stem is None:
            errors.append(f"UNPARSEABLE NAME: {f}")
            continue
        (byfwd if read == "1" else byrev)[stem] = f
    # pair them
    all_stems = set(byfwd) | set(byrev)
    for stem in sorted(all_stems):
        r1 = byfwd.get(stem); r2 = byrev.get(stem)
        if r1 is None or r2 is None:
            errors.append(f"UNPAIRED in {label}: stem '{stem}' R1={r1} R2={r2}")
            continue
        run = flowcell_of(r1)
        sample_id = f"{label}__{stem}"
        runs[run].append((sample_id, r1, r2))
        audit.append((sample_id, label, stem, run, r1, r2))

# fail loudly before writing anything
if errors:
    sys.stderr.write("\n!!! PROBLEMS FOUND — no manifests written:\n")
    for e in errors:
        sys.stderr.write("  " + e + "\n")
    sys.exit(1)

# write per-run manifests
for run, rows in sorted(runs.items()):
    path = os.path.join(OUTDIR, f"manifest_{run}.tsv")
    with open(path, "w") as out:
        out.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for sid, r1, r2 in sorted(rows):
            out.write(f"{sid}\t{r1}\t{r2}\n")

# audit + summary
with open(os.path.join(OUTDIR, "sampleid_to_file_map.tsv"), "w") as out:
    out.write("sample-id\tbatch\toriginal-stem\trun\tR1\tR2\n")
    for row in sorted(audit):
        out.write("\t".join(row) + "\n")

with open(os.path.join(OUTDIR, "run_summary.tsv"), "w") as out:
    out.write("run\tn_samples\tmanifest\n")
    for run, rows in sorted(runs.items()):
        out.write(f"{run}\t{len(rows)}\tmanifest_{run}.tsv\n")

# stdout summary
print(f"Wrote {len(runs)} per-run manifests to {OUTDIR}")
print(f"{'RUN (instr_run_flowcell)':32} {'n_samples':>9}")
print("-"*44)
total=0
for run, rows in sorted(runs.items()):
    print(f"{run:32} {len(rows):>9}")
    total += len(rows)
print("-"*44)
print(f"{'TOTAL':32} {total:>9}")
print(f"\nAudit map: {OUTDIR}/sampleid_to_file_map.tsv")
print(f"Run summary: {OUTDIR}/run_summary.tsv")
