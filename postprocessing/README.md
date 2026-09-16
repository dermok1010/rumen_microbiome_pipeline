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

**Rerun 2026-09-16 against the new data.** Both notebooks originally
hardcoded `os.chdir("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/")`
and read from the **old, pre-rebuild** HPC checkout's directory layout
(`full_combined_results/...`, `corrected_taxonomy_515F_806R/taxonomy_export/...`).
Confirmed by row count at the time: `rumen-core`'s existing
`genus_counts_no_controls_515F_806R.csv` (dated 2026-09-12, 1459 rows)
matched this notebook's own `print` comment ("no_controls: 1459") exactly,
so it was produced by these notebooks against the *old* checkout.
Repointed both at this repo's `../results/rumen_combined_feature-table.tsv`
/ `../results/rumen_combined_exported-taxonomy.tsv`, fixed an
undefined-variable bug in notebook 02 (`prev[keep]` in the feature-summary
cell referenced a `prev` that was never defined after an earlier refactor
to per-host-group prevalence — now computed as overall prevalence across
all samples), and reran. Output goes to `output/` (gitignored) here, **not**
`../results/`.

New counts: 1761 samples total, 1733 after control-filtering (28 flagged,
same regex as before), 238 genus-level features retained at ≥10%
prevalence in at least one host group (Sheep 1126, Beef 164, Dairy 165,
"Other"/unclassified 278 — that host-group split is inherited unchanged
from `02_prepare_model_matrix.ipynb`'s original `assign_species` logic,
not something touched here).

`output/` here has **not** been copied into `rumen-core/data/` — deliberate.
`rumen-core/docs/analysis-scope.md` names the existing
`genus_counts_no_controls_515F_806R.csv`/`CLR_genus_only_515F_806R.csv`
as backing an active paper, and the new sample count (1733 vs 1459) is a
real change to the underlying data, not a like-for-like rerun. Promoting
these files there needs the user's explicit sign-off.

Requires the `rumen_microbiome_pipeline` micromamba env
(`~/envs/rumen_microbiome_pipeline`; registered as Jupyter kernel
`rumen_pp`) — created 2026-09-16 since no existing project env on this VM
had both pandas/numpy and R readr/dplyr/tidyr/readxl.

## `phenotype_merge_scripts/`

Four ad hoc R scripts, originally run interactively on HPC (`setwd()`/
hardcoded `/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/...` paths
throughout — not a packaged pipeline). They map sequenced microbiome
samples to animal phenotype records, one source dataset each. Related to,
but distinct from, the genus-to-CLR step above: this is ID/phenotype
merging (microbiome sample → animal → phenotype), not the count-table
transform.

- **`sheep_order_matching.R`** — FASTQ sample name → tube ID → EID +
  sample date, across 4 sheep studies (CT24, NZAC, Haggot, CTmicro).
  **Not rerun**: two of its four source files (a Grange-lab Excel sheet
  and an old NAS-mounted `Macrogen Sample List Feb 2025.csv`) aren't
  present anywhere on this VM. Its prior output,
  `~/hpc_incoming/sheep_data/sheep_sample_to_EID.csv` (dated 2026-08-31,
  matching the script's own timestamp), was reused as-is instead of
  regenerating it.
- **`sheep_methane_merge.R`** — joins the above to nearest-date methane
  measurements from the PAC pipeline's `PACfile_ani_id.csv`, then aligns
  genus/CLR/ASV data to the matched samples. A more developed VM-adapted
  copy of this script already existed at
  `~/sheep-methane-genomics-microbiome/analysis_work/microbiome_merge_input/02_merge_sheep_methane_vm.R`,
  last run 2026-09-14 against the *old* pre-rebuild genus/CLR data.
  Reran 2026-09-16 with the new 17-flowcell genus/CLR/feature-table data
  in a new sibling directory,
  `analysis_work/microbiome_merge_input_17flowcell/` (old-data output
  left untouched for comparison). Result: identical sample/animal counts
  (1116 samples, 1036 animals) to the old-data run, ASV count shifted
  slightly (38256 vs 38289) from the new pipeline run.
- **`tully_merge.R`** — merges a beef (Tully) genus-count table against
  VFA/pH/methane phenotype data by animal tag. Its original input,
  `tully_data/genus_counts.csv`, no longer exists anywhere on this VM.
  Reconstructed a Tully-only equivalent
  (`phenotype_merge_scripts/output/tully_genus_counts.csv`, gitignored)
  by filtering the new `genus_counts_no_controls_515F_806R.csv` to
  `SampleID` starting with `Tully__` and splitting `ANI_ID`/`SeqID` out
  of the sample name (e.g. `Tully__10732_S4` → ANI_ID `10732`, SeqID
  `S4` — all 126 Tully rows matched this pattern cleanly). Ran the
  merge (`output/run_tully_merge.R`, same logic as `tully_merge.R` with
  `library(dplyr)` moved to the top and the three I/O paths updated) via
  `output/beef_full_clean.csv`: 124 of 124 Tully animals matched a
  phenotype record, 3 dropped as ID-ambiguous in either dataset, 121
  final rows, all with a `CH4__g_d_` value.
- **`charlie_1.R`** (dairy) — **not run**: all 4 of its input files
  (bridge sample-ID sheet, HM24/CRT23/clover22 cow-data sheets) are
  missing from this VM. Skipped per the user's direction 2026-09-16;
  revisit once those files are uploaded to `hpc_incoming`.
