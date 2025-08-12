# Rumen Microbiome Analysis

This project processes 16S rRNA sequencing rumen samples to generate genus-level relative abundances and calculate alpha diversity metrics. Using Qiime package raw fasta files are checked deoinised and matched with SILVA refernce database

Qiime output files were generated and data was also exported to python (jupyter) in which nromalised relative abundances and diversity metrics were calculated

rumen_microbiome_pipeline/
├── create_manifest_file.sh # Generates wide-form manifest for QIIME2
├── qiime_import.slurm # Imports paired-end FASTQs into QIIME2
├── qiime_dad2.slurm # Denoises data and generates ASV table and rep seqs
├── qiime_summarise.slurm # Summarises DADA2 outputs for QC
├── qiime_taxonomy.slurm # Assigns taxonomy using SILVA classifier
├── export.sh # Exports table & taxonomy to TSVs
├── qiime_barplot.slurm # (Optional) Generate taxonomy barplots
├── relative_abundance.ipynb # Python script to calculate relative abundances
