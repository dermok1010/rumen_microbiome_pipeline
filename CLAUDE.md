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

- The genus-to-CLR / "no_controls" filtering post-processing step used by
  `rumen-core` and the sheep-methane Paper 4 merge is **not in this
  repository** -- `collapse_genus` produces raw genus counts only.
  Whatever produces the `_no_controls`/`CLR_genus_only` files wasn't
  found in the source snapshot this rebuild is based on. Track it down
  before assuming this repo alone reproduces those downstream files.
- **Nextflow (26.04.6)/Java 17/Docker are installed on this VM** (2026-09-14)
  and `main.nf`/`nextflow.config` are syntax-validated -- both entries
  parse and complete with zero tasks against no data. Requires
  `NXF_SYNTAX_PARSER=v1` (newer Nextflow's strict parser rejects this
  pipeline's classic multi-`workflow` + `-entry` pattern); `run/run_nextflow.slurm`
  sets it, set it yourself for manual invocations. Not yet validated: an
  actual run against real/fixture data, or running (not just building)
  the Docker image.
- The classifier's SILVA reference assets and trained classifier binary
  aren't tracked (large binaries, by design) -- see
  `classifier/assets/README.md`.
- This repo is not yet wired into the Claude-worktree system the other
  three projects use (`agent-worktrees/claude/...`) -- it's a plain
  checkout for now, same pattern as `sheep-methane-genomics-microbiome`.
