# rumen_microbiome_pipeline

16S rRNA amplicon processing (QIIME2 + DADA2, Nextflow DSL2) that turns
raw rumen microbiome FASTQs (sheep, beef, dairy) into ASV/genus count
tables. This is the *upstream* pipeline -- it is not `rumen-core`.
`rumen-core`'s own `CLAUDE.md` describes itself as downstream analysis
"using genus-level count/CLR data produced upstream on HPC"; this repo is
that upstream producer. See `README.md` for pipeline stages and the
classifier/primer-correction history.

See also `~/.claude/CLAUDE.md` for VM-wide conventions.

## History -- rebuilt 2026-09-14

This repo was substantially rebuilt from a source snapshot
(`rumen_microbiome_pipeline_source_2026-09-14.tar.gz`, from the HPC path
`/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline`,
~123GB full size on HPC) rather than from what was previously on GitHub --
GitHub's prior version was a single-batch `main.nf` that predated the
per-flowcell batch-aware rewrite and the 515F/806R primer correction
entirely. The full audit that preceded this rebuild is preserved at
`~/staging/rumen_microbiome_pipeline_audit/AUDIT_REPORT.md` on this VM
(not in git -- it documents the archaeology, not the pipeline itself).

The **HPC copy is being kept as-is, untouched, as a backup** of the state
before this rebuild. This repo and the HPC's existing checkout will
diverge until the user deliberately re-deploys the rebuilt version there.

## Directory layout

- `main.nf` / `nextflow.config` -- the pipeline. Two entrypoints:
  `import_only`, `full_run`. `nextflow.config`'s `standard` profile
  (default) is local-executor, no container, for VM/fixture use; `hpc`
  is Slurm + Apptainer for production.
- `make_manifests.py` -- builds per-flowcell manifests from raw FASTQs.
  Supports `--batch`/`--flowcell`/`--sample-id-file` to generate a
  deliberate subset rather than requiring a full manifest edit.
- `classifier/build/` -- the *active* 515F/806R classifier-building
  scripts (extract reference reads, train). `classifier/legacy_341F_806R/`
  is the superseded wrong-primer version, kept for provenance only, not
  used by anything. `classifier/assets/` (gitignored) is where the large
  binary `.qza` files actually live once built.
- `run/` -- HPC Slurm wrappers (`run_manifests.slurm`, `run_nextflow.slurm
  <entry>`). These are inherently HPC-specific (module loads, Slurm
  directives) -- that's correct, don't try to make them VM-portable.
- `legacy_manual_pipeline/` -- pre-Nextflow standalone scripts. Already
  archived by the pipeline's original author (commit message: "Refactor
  repo: promote Nextflow pipeline and archive legacy manual workflow");
  carried forward unchanged, reference value only.

## Data & results

- `data/`, `outputs/`, `results/`, `manifest_files/` are VM-only and
  gitignored -- never commit them, never assume they exist in a fresh
  checkout.
- **No fixture data exists on the VM yet.** `data/raw/` is where a small
  representative FASTQ subset should eventually go for local/VM testing.
  Extracting one is separate, deliberate future work -- don't invent one
  ad hoc; choose it based on the flowcell/species breakdown in the audit
  report.

## Known gaps (carried forward from the audit, not yet resolved)

- **Rerun 2026-09-16 against the new 17-flowcell data.** The genus-to-CLR
  / "no_controls" filtering step (`collapse_genus` in `main.nf` still
  produces raw genus counts only) is handled by
  `postprocessing/01_build_microbiome_matrices.ipynb` and
  `02_prepare_model_matrix.ipynb` -- both repointed from the old
  pre-rebuild HPC paths to this repo's `results/` layout and rerun.
  Fixed an undefined-variable bug in notebook 02 (`prev[keep]` in the
  feature-summary cell referenced a `prev` that was never defined after
  an earlier refactor to per-host-group prevalence; now computed as
  overall prevalence across all samples). Output lands in
  `postprocessing/output/` (gitignored, not `results/`) -- **this has
  NOT been copied into `rumen-core/data/`**, deliberately: the new run
  covers 1761 samples (1733 after control-filtering) vs the old
  1487/1459, a real change in the underlying data, and
  `rumen-core/docs/analysis-scope.md` names the existing
  `genus_counts_no_controls_515F_806R.csv`/`CLR_genus_only_515F_806R.csv`
  as backing an active paper. **The user has confirmed (2026-09-16) they
  intend to eventually promote these new files into `rumen-core/data/`
  and the sheep-methane Paper 4 merge, replacing the older versions --
  just not in this pass.** Don't do it proactively; wait for them to ask.
  See `rumen-core/docs/analysis-scope.md`'s "Pending data refresh" section
  for the specifics (notably: the Sheep/Beef/Dairy host-group counts came
  out identical to `rumen-core`'s hardcoded reproducibility-guard values,
  so that side is expected to be low-risk; the sheep-methane side already
  reran cleanly with matching sample/animal counts -- see that repo's
  `analysis_work/microbiome_merge_input_17flowcell/`). Requires the
  `rumen_microbiome_pipeline`
  micromamba env (`~/envs/rumen_microbiome_pipeline` -- python
  pandas/numpy/nbclient + R readr/dplyr/tidyr/readxl; registered as
  Jupyter kernel `rumen_pp`), created 2026-09-16 as this project's own
  env since no existing project env covered both notebook and R-script
  dependencies without borrowing across projects.
- **Nextflow (26.04.6)/Java 17/Docker are installed on this VM.** Requires
  `NXF_SYNTAX_PARSER=v1` (newer Nextflow's strict parser rejects this
  pipeline's classic multi-`workflow` + `-entry` pattern); `run/run_nextflow.slurm`
  sets it, set it yourself for manual invocations. There's no local qiime2
  install, so the `standard` profile can't run here -- use the
  `docker_local` profile instead (2026-09-15), which runs against
  `ghcr.io/dermok1010/rumen-microbiome-pipeline:latest`, the same image
  `build_sif.slurm` pulls for the HPC `.sif`. Verified end-to-end
  (`import_only` against a synthetic FASTQ pair, not just a config-parse
  check) -- real production data still hasn't been run through `docker_local`
  specifically. The `hpc` profile has: a 17-flowcell production run
  completed on HPC 2026-09-15/16 and its output (`outputs/`, `results/`)
  was pulled back to this VM -- see the "Known gaps" genus-to-CLR entry
  below for what still needs to happen to it.
- The classifier's SILVA reference assets and trained classifier binary
  aren't tracked (large binaries, by design) -- see
  `classifier/assets/README.md`. On the HPC checkout, the trained
  515F/806R classifier was copied in from the pre-rebuild
  `Paper_2/rumen_microbiome_pipeline` checkout (2026-09-15) rather than
  rebuilt -- same file, no need to redo the SILVA extract+train.
- This repo is not yet wired into the Claude-worktree system the other
  three projects use (`agent-worktrees/claude/...`) -- it's a plain
  checkout for now, same pattern as `sheep-methane-genomics-microbiome`.
