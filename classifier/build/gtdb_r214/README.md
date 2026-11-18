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

## Full reclassification for rumen-core (VM, 2026-11-18)

rumen-core's network paper moved from SILVA 138 to GTDB r214, so every ASV of the combined run is
now classified, not only the comparison animals:

    bash classifier/build/gtdb_r214/05_classify_all_asvs.sh   # export, filter, classify, then 06
    python classifier/build/gtdb_r214/06_full_genus_table.py  # genus tables, reconciliation, crosswalk

Output: `outputs/gtdb_r214_full/` (gitignored). `05` sets `OUT` before sourcing `common.sh`, which
now defaults `OUT` to the old folder only when it is unset, so `03`/`04` behave as before. About 6
minutes on 4 cores.

- **Inputs.** The same ASV table and SILVA taxonomy as `postprocessing/01_build_microbiome_matrices.ipynb`
  (`results/rumen_combined_feature-table.tsv` and `rumen_combined_exported-taxonomy.tsv`, 1 Oct 2026
  run), which made `postprocessing/output/genus_counts_no_controls_515F_806R.csv`. Rebuilding that
  SILVA table from these files with the notebook's rules reproduces it exactly (max cell difference 0,
  1,919 samples), so the GTDB table has the same samples and per-sample totals.
- **Rep seqs.** Only the 15 Sep rep seqs are on the VM. 6,774 of the 90,464 ASVs with reads in
  non-control samples have no sequence (Tully ASVs new in the 1 Oct rerun): 0.10% of all
  non-control reads (up to 20% in one low-depth sample outside rumen-core), and 0.0017% of the
  rumen-core phenotyped subset's reads (max 0.34% in one sample, beef only). They are pooled as
  `Unassigned`.
- **Genus handling** mirrors notebook 01: the `g__` rank without its prefix; no genus rank, an
  empty one, or no taxonomy -> one `Unassigned` column. GTDB suffixes (`_A`, `_B`) stay distinct
  genera. Controls are removed with the notebook's two regular expressions.
- **Result.** 1,500 GTDB genera (SILVA: 924 labels, some of them `uncultured`/`Unknown_Family`).
  Reads assigned to a genus in the rumen-core subset: 79.4% GTDB vs 95.1% SILVA label (93.5% excluding
  SILVA's `uncultured`/`Unknown_Family`); the gap is GTDB genome-only genera with no 16S reference.
- **Crosswalk** `silva_gtdb_genus_crosswalk.csv`: for each SILVA genus, the GTDB genera its ASVs'
  reads go to, with read shares (rumen-core subset and all non-control samples) and ASV counts.
  E.g. *Methanobrevibacter* -> *Methanobrevibacter_A* 77%, *Methanobrevibacter* 22%;
  *Shuttleworthia* -> *Shuttleworthia* 40%, Unassigned 33%, *Weimeria* 16%; *Kandleria* ->
  *Catenibacterium* 99.9%; Christensenellaceae R-7 group -> RUG740 63%, Unassigned 18%.
