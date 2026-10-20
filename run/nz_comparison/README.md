# Stand-alone 16S run: NZ_comparison (16S vs RE-RRS samples)

Processes BSI's 16S run of the `comparison_16s_vs_rerrs_2024` samples
(92 NZ Flock 2638 sheep, 92 Irish INZAC sheep, 164 Irish cattle, plus
controls/replicate libraries; ~361 libraries) through this pipeline,
replacing the one-off 2026-09 QIIME2 table that
`micro_prediction_NZ/data/re_rrs/comparison_16s_vs_rerrs_2024/16S_table/table.qza`
holds. It is kept out of the combined V4 run because it is a separate
sequencing run and probably a different amplicon: that table's QIIME2
provenance shows `trim_left 17/21`, `trunc_len 280/260`, which matches
341F/805R (V3-V4, 2x300) rather than this pipeline's 515F/806R (V4).
Step 1 confirms this before anything is submitted.

Everything runs in its own directory (`NZ_RUN_DIR`) so the combined run's
`outputs/`, `results/`, `work/` and `-resume` history are untouched.

## Steps (on the HPC)

```bash
PIPE_DIR=/path/to/rumen_microbiome_pipeline        # HPC checkout
NZ_RUN_DIR=/path/to/nz_comparison_16s_run          # new, empty
cd "$PIPE_DIR" && git pull                          # commit with run/nz_comparison/
```

1. **Check the reads** (login node, ~1 min, no Slurm):
   `python3 run/nz_comparison/check_reads.py`
   (default folder `/data/Genetics/primary/R1240_microbiome/NZ_comparison`;
   pass another path if the FASTQs live elsewhere, e.g. under
   `/data/BioScience/primary/R2002_methanepredict/NZ_method_comparison`).
   Note the primer fractions, read length and quality-by-position.
2. **Fill the params**: `mkdir -p "$NZ_RUN_DIR" && cp run/nz_comparison/nz_comparison.params.yaml "$NZ_RUN_DIR/"`,
   replace `PIPE_DIR`/`NZ_RUN_DIR`, then:
   - 341F/805R reads -> keep `trim_left 17/21` and the 341F-805R classifier;
   - 515F/806R reads -> `trim_left 19/20` and `silva-138-99-515F-806R-classifier.qza`;
   - set `trunc_len_f/r` just before median quality falls below ~25, keeping
     enough overlap (see the comment in the params file).
   If the data folder is not `R1240_microbiome/NZ_comparison`, change
   `OPT_IN_BATCHES` in `make_manifests.py` (commit it) before step 3.
3. **Manifests** (login node or `run_manifests.slurm`-style env; stdlib only):
   `python3 make_manifests.py manifest_files/nz_comparison --batch NZcmp`
   Check `manifest_files/nz_comparison/run_summary.tsv`: expect one or a
   few flowcells and ~361 sample-ids (`NZcmp__<stem>`).
4. **Classifier** (only for 341F/805R reads, and only if
   `classifier/assets/silva-138-99-341F-805R-classifier.qza` doesn't exist):
   `cd classifier/assets && sbatch ../../run/nz_comparison/build_classifier_341F_805R.slurm`
5. **Pipeline**, after the classifier job finishes:
   `cd "$NZ_RUN_DIR" && PIPE_DIR="$PIPE_DIR" sbatch "$PIPE_DIR/run/nz_comparison/run_nz_comparison.slurm"`
6. **Bring results back** to the VM (current VM IP via `curl -s ifconfig.me`
   on the VM), into `~/hpc_incoming/nz_comparison_16s/`: `results/`,
   `outputs/*.qzv`, `outputs/per_batch/*/*_denoising_stats.qzv`,
   `run_commit.txt`, `run_params.yaml`.

## Rough cost (Background queue)

| Job | Resources | Expected wall time |
|---|---|---|
| Classifier build (step 4) | 4 CPU, 30 GB | ~3-6 h (extract + naive Bayes fit on SILVA 138.99) |
| Nextflow driver (step 5) | 2 CPU, 8 GB | for the whole run |
| import + DADA2 per flowcell | 8 CPU, 32 GB | ~2-6 h for ~361 libraries |
| Taxonomy, genus collapse, export | 8 CPU, 24 GB | ~1 h |

Estimates are by analogy with the combined run, not measured for this data.

## After the run

Sample-ids are `NZcmp__<stem>`; the old table used the bare `LabSampleID`
(e.g. `935721`, `935721a`). Rebuild
`micro_prediction_NZ/build_16s_re_rrs_sample_crosswalk.py`'s join on the
new ids (strip the prefix and any `_S<N>` tag) before using the table.
