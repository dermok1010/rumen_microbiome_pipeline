# Classifier reference assets (not tracked in git)

This directory holds the large binary QIIME2 artifacts the classifier
build/use steps need. All `.qza` files here are gitignored (SILVA
reference release is large, and the trained classifier itself is a
multi-hundred-MB binary) -- fetch/build them once per environment:

- `silva-138-99-seqs.qza`, `silva-138-99-tax.qza` -- SILVA 138.99
  reference sequences/taxonomy (download once from the QIIME2 data
  resources page, or copy from the HPC's existing copy).
- `silva-138-99-515F-806R.qza` -- built by `classifier/build/extract_reads_515F_806R.slurm`.
- `silva-138-99-515F-806R-classifier.qza` -- built by
  `classifier/build/train_classifier_515F_806R.slurm`. This is
  `nextflow.config`'s default `params.classifier`.

Submit both build scripts from inside this directory (`sbatch
extract_reads_515F_806R.slurm` etc., run from `classifier/assets/`, not
`classifier/build/`).
