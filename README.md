# rumen_microbiome_pipeline

16S rRNA amplicon processing (QIIME2 + DADA2) of rumen microbiome samples
(sheep, beef, dairy). Import per-sequencing-flowcell -> DADA2 denoise
(separate error model per flowcell) -> merge -> taxonomy classification
-> genus-level collapse -> export.

This is the upstream pipeline that produces the genus-count and ASV
feature tables consumed by the `rumen-core` project and the sheep-methane
paper's microbiome-extension work.

## Pipeline

```
import_reads (per flowcell) -> dada2_denoise (per flowcell)
  -> merge_batches -> summarise_outputs
                    -> assign_taxonomy -> collapse_genus
                    -> export_outputs
```

Two Nextflow entrypoints in `main.nf`:
- `import_only` -- import + demux-summarize every flowcell, for eyeballing
  read quality before committing to a full denoise.
- `full_run` -- the full pipeline above.

Batching is by sequencing flowcell (read from each FASTQ's own header via
`make_manifests.py`), not by species -- DADA2 needs a separate error model
per flowcell. This is also the origin of the `<BATCH>__<stem>` sample-ID
convention (e.g. `EN00011687__Beef_2023_31`, `CTmicro__00794`) used
throughout the other three projects on this VM.

## Classifier

The pipeline classifies against a SILVA 138.99 naive-Bayes classifier
trained on **515F/806R**-extracted reference reads (`classifier/build/`),
matching the actual wet-lab PCR primers. An earlier version of this
pipeline shipped with a classifier trained on the wrong primer region
(341F/806R) -- see `classifier/legacy_341F_806R/README.md`. That path is
not supported; if the classifier ever needs rebuilding (new SILVA
release, etc.), extend `classifier/build/`, don't fork a second primer
path.

The trained classifier itself is a large binary and is gitignored --
build it once per environment following `classifier/assets/README.md`.

## Running it

```bash
# VM / local / fixture testing (default profile, local executor, no container):
python3 make_manifests.py manifest_files/combined_run --data-root /path/to/fixture/data
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -entry full_run

# A deliberate subset (by batch, flowcell, or sample list) instead of everything:
python3 make_manifests.py manifest_files/subset --batch EN00011687 --flowcell M07726_229

# HPC production run (Slurm + Apptainer):
sbatch run/run_manifests.slurm
sbatch run/run_nextflow.slurm full_run     # or: import_only
```

`NXF_SYNTAX_PARSER=v1` is required: newer Nextflow defaults to a stricter
parser that rejects this pipeline's classic multiple-`workflow`-block +
`-entry` pattern ("the -entry option is not supported with the strict
parser"). `run/run_nextflow.slurm` already sets it; set it yourself for
any manual/interactive `nextflow run` invocation.

All paths (`raw_data_root`, `manifest_dir`, `classifier`, `output_dir`,
`results_dir`, Apptainer bind-mounts) are `params.*` in `nextflow.config`
with VM-friendly defaults -- override per-invocation with `--param value`
or a `-params-file` for the HPC's actual paths. Nothing in the active
pipeline hardcodes `/home/dermot.kelly/...` or `/data/Genetics/...` as the
only option; those remain the HPC's own override values, not baked-in
defaults.

## Intended workflow

Develop and test on the VM against a small fixture -> commit the exact
version to GitHub -> `git pull` the same commit on the HPC -> run against
the full production dataset with `-profile hpc` and the HPC's own
params/bind-mounts -> record the commit hash, config, container digest,
and exact manifest/command alongside the results.

## Known gaps

- **No fixture data exists yet.** `data/raw/` (gitignored) is where a
  small representative FASTQ subset should eventually live for VM-side
  testing; it hasn't been extracted from the HPC yet.
- **The genus-to-CLR / "no_controls" filtering step used downstream
  (`rumen-core`, the sheep-methane Paper 4 merge) is not part of this
  repository.** This pipeline's `collapse_genus` process produces raw
  genus counts only. Whatever script filters out control samples and
  computes the CLR transform lives elsewhere and wasn't found in the
  source snapshot this rebuild was based on -- track it down before
  relying on this repo alone to reproduce those downstream files.
- **Nextflow (26.04.6), Java 17, and Docker are now installed on this VM**
  and `main.nf`/`nextflow.config` have been syntax-checked (`nextflow
  config`, `nextflow run -entry import_only`/`full_run` with no data --
  both parse and complete with zero tasks, as expected with no fixture
  present yet). Not yet validated: an actual end-to-end run against real
  or fixture FASTQ data, and the built Docker image hasn't been run
  (only built) -- confirm both once a fixture exists.
- The classifier's own reference assets (SILVA release, trained
  classifier binary) aren't in this repo by design (large binaries) --
  see `classifier/assets/README.md` for how to (re)build them.
