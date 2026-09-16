# Post-processing: genus-to-CLR and phenotype merging

Staged here 16 Sept 2026 from `~/hpc_incoming/` (pushed from HPC). These
recover the step flagged as missing in this repo's `CLAUDE.md` "Known gaps"
section since the 2026-09-14 rebuild: `main.nf`'s `collapse_genus` process
produces raw genus counts only, and the control-filtering/CLR transform that
`rumen-core` and the sheep-methane Paper 4 merge actually consume wasn't in
the source snapshot the rebuild was based on.

## `01_build_microbiome_matrices.ipynb` → `02_prepare_model_matrix.ipynb`

Run in that order. Together they take the QIIME2 export (feature table +
taxonomy) and produce the two files `rumen-core` names as its "current
processed inputs" (`docs/analysis-scope.md`):

1. **`01_build_microbiome_matrices.ipynb`** — collapses the ASV feature
   table to genus level, flags and strips control/blank samples by regex
   (`__(?:ControlN_ive|ControlNive|ControlP_ive|Minus|Plus|N|P)_?\d*$`),
   writes `genus_counts_with_controls_515F_806R.csv` and
   `genus_counts_no_controls_515F_806R.csv`, plus `taxonomy_lookup_515F_806R.csv`.
2. **`02_prepare_model_matrix.ipynb`** — loads the no-controls genus table,
   derives a Sheep/Beef/Dairy host group from `SampleID`, keeps taxa at
   ≥10% prevalence within at least one host group, CLR-transforms
   (pseudocount `1e-6`), and writes `CLR_genus_only_515F_806R.csv`,
   `raw_filtered_genus_515F_806R.csv`, `feature_summary_515F_806R.csv`.

**Path mismatch, not yet fixed:** both notebooks hardcode
`os.chdir("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/")`
and read from `full_combined_results/rumen_combined_feature-table.tsv` and
`corrected_taxonomy_515F_806R/taxonomy_export/taxonomy.tsv` — the **old,
pre-rebuild** HPC checkout's directory layout (`Paper_2/...`), not the
rebuilt pipeline's actual `results/rumen_combined_feature-table.tsv` /
`results/rumen_combined_exported-taxonomy.tsv` that landed in `../outputs/`
and `../results/` from the 17-flowcell run on 2026-09-15/16. Confirmed by
row count: `rumen-core`'s existing `genus_counts_no_controls_515F_806R.csv`
(dated 2026-09-12) has exactly 1459 rows, matching this notebook's own
`print` comment ("no_controls: 1459") — so these notebooks produced that
existing copy from the *old* checkout, not the new run. Before rerunning
against `../results/`, update the `FEATURE_FILE`/`TAX_FILE`/`os.chdir` paths
in both notebooks, and expect the sample count to change (17 flowcells here
vs whatever the old checkout covered).

**Not yet done:** the notebooks have not been rerun against the new
`../results/` output, and `rumen-core/data/{genus_counts_no_controls,
CLR_genus_only}_515F_806R.csv` have not been refreshed. `rumen-core`'s own
`docs/analysis-scope.md` states these files back an active paper
("Reproducibility principle: the first VM run must reproduce the current
downstream analysis as received from HPC") — don't regenerate/overwrite
`rumen-core`'s copies without the user's sign-off, since that changes the
scientific data behind that manuscript's current state.

## `phenotype_merge_scripts/`

Four ad hoc R scripts, run interactively on HPC (`setwd()`/hardcoded
`/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/...` paths throughout —
not a packaged pipeline). They map sequenced microbiome samples to animal
phenotype records, one source dataset each:

- `sheep_order_matching.R` — FASTQ sample name → tube ID → EID + sample
  date, across 4 sheep studies (CT24, NZAC, Haggot, CTmicro).
- `sheep_methane_merge.R` — joins the above to nearest-date methane
  measurements from the PAC pipeline's `PACfile_ani_id.csv`.
- `tully_merge.R` — merges a beef (Tully) genus-count table against
  VFA/pH/methane phenotype data by animal tag.
- `charlie_1.R` — cross-checks which dairy cow-data sheets (HM24/CRT23/
  clover22) resolve a bridge sample-ID sheet's farm names.

Related to, but distinct from, the genus-to-CLR step above: this is
ID/phenotype merging (microbiome sample → animal → phenotype), not the
count-table transform. Kept for provenance; not adapted to run here yet.
