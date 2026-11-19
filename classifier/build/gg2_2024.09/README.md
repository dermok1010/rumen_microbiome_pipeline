# Greengenes2 2024.09 classification of all ASVs

**Why.** PhD Chapter 5 (`~/methane-signature-replication`, private repo
`dermok1010/methane-signature-replication`) projects our cohorts onto the published NMF
"Ruminosignatures" of Vourlaki et al. 2026 (ISME J, doi 10.1093/ismejo/wrag185). Their W matrices are
indexed by Greengenes2 genus lineages. With GTDB r214 names, UBA2810 (about 60% of the headline
UBA2 signature) and *Kandleria* had no match, and several GTDB/GG2 name splits lowered fidelity.
Their methods state taxonomy came from "a customized GreenGenes2 (release 2024.09)" trained for each
dataset's primer set; their Irish data are V4. This folder replicates that route on our V4 ASVs.

This does **not** replace the pipeline's SILVA classifier or the GTDB r214 tables; it is a third
reference used for the Ch5 replication only.

## Steps (VM, pipeline Docker image, run from the repo root)

    bash classifier/build/gg2_2024.09/01_get_reference.sh      # GG2 2024.09 backbone -> ~/refs/greengenes2_2024.09 (+ MD5SUMS, SOURCE.txt)
    bash classifier/build/gg2_2024.09/02_train_classifier.sh   # naive Bayes on backbone V4 (515F/806R) + backbone taxonomy
    bash classifier/build/gg2_2024.09/03_classify_all_asvs.sh  # classify-sklearn all ASVs, then 04
    ~/envs/rumen_microbiome_pipeline/bin/python classifier/build/gg2_2024.09/04_genus_table.py

`03` reuses `outputs/gtdb_r214_full/all_rep_seqs.qza` (the 83,690 ASVs with reads in a non-control
sample that have a rep sequence), so `gtdb_r214/05_classify_all_asvs.sh` (branch
`gtdb-r214-full-reclassification`, PR #1) must have been run first.

## Choices

- **Own classifier, not GG2's pre-trained `2024.09.backbone.v4.nb.qza`.** That file is a pickled
  scikit-learn object, i.e. code; training from the backbone FASTA/taxonomy (plain data) avoids
  executing a download and guarantees the scikit-learn version matches the image's. Default
  q2-feature-classifier settings; no dereplication (class priors as in GG2's own model);
  `--p-classify--chunk-size 2000` only to bound memory (28,054 taxonomy strings).
- **Confidence 0.7** (classify-sklearn default), same as the GTDB r214 run.
- **Genus columns are full lineages** (`Bacteria;Pseudomonadota;...;Succinivibrionaceae;UBA2810`),
  prefixes stripped, GG2 numeric suffixes kept: identical to the W row names.

## Where things go

| What | Where | In git |
|---|---|---|
| Scripts, this README | `classifier/build/gg2_2024.09/` | yes |
| GG2 release files, MD5SUMS, SOURCE.txt | `~/refs/greengenes2_2024.09/` | no |
| Trained classifier | `classifier/assets/gg2_2024.09/` | no |
| ASV taxonomy, genus table, summary | `outputs/gg2_2024.09_full/` | no |
