# GTDB r214 classifier for the 16S vs RE-RRS comparison

**Why.** BSI's RE-RRS genus counts use GTDB r214 (sheep; Tully beef RE-RRS is r220), while this
pipeline classifies against SILVA 138.99. The two taxonomies only share ~half the reads by genus
name, so genus-by-genus comparisons between the technologies need the 16S ASVs classified on
GTDB too. Decision (Dermot, 2026-10-27): GTDB **r214 for both species**.

This does **not** replace the pipeline's SILVA classifier (`nextflow.config`'s
`params.classifier`); it is a second reference used for the comparison only. Same primers
(515F/806R) and the same naive-Bayes method as `../extract_reads_515F_806R.slurm` and
`../train_classifier_515F_806R.slurm`, so only the reference changes.

## Where things go

| What | Where | In git |
|---|---|---|
| Scripts, this README | `classifier/build/gtdb_r214/` | yes |
| GTDB reference, V4 extract, trained classifier | `classifier/assets/gtdb_r214/` | no |
| ASV taxonomy and genus tables | `outputs/gtdb_r214_reclassification/` | no |

## Steps (VM, via the pipeline's own Docker image; run from the repo root)

    bash classifier/build/gtdb_r214/01_get_reference.sh    # download GTDB 214.1 SSU, species reps
    bash classifier/build/gtdb_r214/02_extract_train.sh    # 515F/806R extract, dereplicate, train
    bash classifier/build/gtdb_r214/03_classify_asvs.sh    # classify the comparison ASVs
    python classifier/build/gtdb_r214/04_genus_tables.py   # genus count tables, sheep and beef

`common.sh` holds the paths and the Docker wrapper (image
`ghcr.io/dermok1010/rumen-microbiome-pipeline:latest`, QIIME2 amplicon 2024.2 with RESCRIPt).

## Choices

- **GTDB SpeciesReps, both domains** (RESCRIPt's default): one 16S sequence per GTDB species.
  Lighter to train than GTDB's "All" set and enough to place a V4 read at genus level.
- **Only the ASVs present in the comparison animals** (505 sheep + 116 Tully beef from
  `micro_prediction_NZ/re_rrs_comparison/prep_inputs/`) are classified: 28,992 ASVs.
- **Representative sequences** come from `outputs/rumen_combined_rep_seqs.qza` (15 Sep 2026 run).
  The ASV table is from the 1 Oct rerun, whose rep seqs stayed on the HPC; DADA2 feature IDs are
  sequence hashes, so the Sep file covers all but 60 of the 28,992 ASVs (0.00% of sheep reads,
  0.02% of beef reads). Those 60 are left unclassified.
- GTDB genera with no 16S sequence in GTDB (many CAG-/RUG- genome genera) cannot be assigned by
  any 16S classifier; their reads stop at family or above.

## First run (VM, 2026-10-27)

Whole chain about 8 minutes on 4 cores; training peaked well inside 16 GB. GTDB 214.1
species reps -> 515F/806R extract -> dereplicated -> 19 MB classifier.

| | Genera | Reads to genus | 16S genera also in RE-RRS | Reads in shared genera (16S / RE-RRS) |
|---|---|---|---|---|
| Sheep (RE-RRS on r214) | 789 | 77.4% | 783 | 77.4% / 84.7% |
| Beef (RE-RRS on r220) | 751 | 81.0% | 708 | 76.6% / 76.6% |

With 16S on SILVA, name matching covered only ~51% of sheep 16S reads and ~41% of RE-RRS reads.
16S percentages are of all reads, including the ~20% not classified to genus.
