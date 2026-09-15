// main.nf  — batch-aware (per-flowcell DADA2, then merge)
nextflow.enable.dsl=2

// --------------------------------------------------
// IMPORT READS  (per flowcell batch)
// --------------------------------------------------
process import_reads {

    tag "${batch} - import"

    publishDir "${params.output_dir}/per_batch/${batch}", mode: 'copy'

    cpus 2
    memory '8 GB'

    input:
    tuple val(batch), path(manifest_file)

    output:
    tuple val(batch), path("${batch}_paired_end_demux.qza"), emit: demux_qza
    path("${batch}_paired_end_demux.qzv")

    script:
    """
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path $manifest_file \
      --output-path ${batch}_paired_end_demux.qza \
      --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
      --i-data ${batch}_paired_end_demux.qza \
      --o-visualization ${batch}_paired_end_demux.qzv
    """
}

// --------------------------------------------------
// DADA2  (per flowcell batch — separate error model each)
// --------------------------------------------------
process dada2_denoise {

    tag "${batch} - dada2"

    publishDir "${params.output_dir}/per_batch/${batch}", mode: 'copy'

    cpus 8
    memory '24 GB'

    input:
    tuple val(batch), path(demux_qza)

    output:
    path("${batch}_table.qza"),     emit: table
    path("${batch}_rep_seqs.qza"),  emit: repseqs
    path("${batch}_denoising_stats.qza"), emit: stats
    path("${batch}_denoising_stats.qzv")

    script:
    """
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $demux_qza \
      --p-trim-left-f ${params.trim_left_f} \
      --p-trim-left-r ${params.trim_left_r} \
      --p-trunc-len-f ${params.trunc_len_f} \
      --p-trunc-len-r ${params.trunc_len_r} \
      --p-n-threads ${task.cpus} \
      --o-table ${batch}_table.qza \
      --o-representative-sequences ${batch}_rep_seqs.qza \
      --o-denoising-stats ${batch}_denoising_stats.qza \
      --verbose

    qiime metadata tabulate \
      --m-input-file ${batch}_denoising_stats.qza \
      --o-visualization ${batch}_denoising_stats.qzv
    """
}

// --------------------------------------------------
// MERGE per-batch tables and rep-seqs into one
// --------------------------------------------------
process merge_batches {

    tag "${params.prefix} - merge"

    publishDir params.output_dir, mode: 'copy'

    cpus 2
    memory '16 GB'

    input:
    path(tables)
    path(repseqs)

    output:
    path("${params.prefix}_table.qza"),    emit: table
    path("${params.prefix}_rep_seqs.qza"), emit: repseqs

    script:
    // build --i-tables / --i-data args from all staged files
    def table_args   = tables.collect  { "--i-tables $it" }.join(' ')
    def repseq_args  = repseqs.collect { "--i-data $it"  }.join(' ')
    """
    qiime feature-table merge \
      ${table_args} \
      --o-merged-table ${params.prefix}_table.qza

    qiime feature-table merge-seqs \
      ${repseq_args} \
      --o-merged-data ${params.prefix}_rep_seqs.qza
    """
}

// --------------------------------------------------
// SUMMARIES  (on merged artifacts)
// --------------------------------------------------
process summarise_outputs {

    tag "${params.prefix} - summary"

    publishDir params.output_dir, mode: 'copy'

    cpus 1
    memory '4 GB'

    input:
    path table
    path repseqs

    output:
    path("${params.prefix}_table.qzv"),    emit: summary_table
    path("${params.prefix}_rep_seqs.qzv"), emit: summary_repseqs

    script:
    """
    qiime feature-table summarize \
      --i-table $table \
      --o-visualization ${params.prefix}_table.qzv

    qiime feature-table tabulate-seqs \
      --i-data $repseqs \
      --o-visualization ${params.prefix}_rep_seqs.qzv
    """
}

// --------------------------------------------------
// TAXONOMY  (once, on merged rep-seqs)
//
// params.classifier must be trained on the SAME primer region the FASTQs
// were amplified with. This pipeline was previously misconfigured to
// classify against a 341F/806R-trained classifier while the actual PCR
// used 515F/806R primers (see classifier/legacy_341F_806R/README.md) --
// the default in nextflow.config now points at the correct 515F/806R
// classifier. There is no supported "other primer set" mode; if the
// wet-lab primers ever change, build a new classifier with
// classifier/build/ and update params.classifier, don't fork this process.
// --------------------------------------------------
process assign_taxonomy {

    tag "${params.prefix} - taxonomy"

    publishDir params.output_dir, mode: 'copy'

    // The one genuinely serial step: runs once on the full merged rep-seqs
    // after every flowcell finishes, so unlike import_reads/dada2_denoise
    // it can't be sped up by running more flowcells in parallel. classify-
    // sklearn's --p-n-jobs parallelizes over this, so more cores here
    // buys real wall-clock time on the critical path. Memory bumped
    // alongside it for the extra parallel worker processes.
    cpus 8
    memory '24 GB'

    input:
    path repseqs
    path classifier_file

    output:
    path("${params.prefix}_taxonomy.qza"), emit: taxonomy
    path("${params.prefix}_taxonomy.qzv")

    script:
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier $classifier_file \
      --i-reads $repseqs \
      --o-classification ${params.prefix}_taxonomy.qza \
      --p-n-jobs ${task.cpus}

    qiime metadata tabulate \
      --m-input-file ${params.prefix}_taxonomy.qza \
      --o-visualization ${params.prefix}_taxonomy.qzv
    """
}

// --------------------------------------------------
// GENUS COLLAPSE  (ASV table -> genus-level count table)
//
// Folds what used to be a manual rerun_taxonomy_*.slurm bolt-on step into
// the pipeline itself, so a fresh run produces genus-level counts
// directly rather than needing a second manual pass. NOTE: this produces
// raw genus counts only. The further "no_controls" filtering and CLR
// transform used downstream (rumen-core, Paper 4) are NOT part of this
// pipeline's tracked source -- that post-processing script lives
// elsewhere and was not found in this repo; see README "Known gaps".
// --------------------------------------------------
process collapse_genus {

    tag "${params.prefix} - genus collapse"

    publishDir params.results_dir, mode: 'copy'

    cpus 1
    memory '8 GB'

    input:
    path table
    path taxonomy

    output:
    path("${params.prefix}_genus_table.qza"), emit: genus_table
    path("${params.prefix}_genus_counts.tsv")

    script:
    """
    qiime taxa collapse \
      --i-table $table \
      --i-taxonomy $taxonomy \
      --p-level 6 \
      --o-collapsed-table ${params.prefix}_genus_table.qza

    mkdir -p genus_export
    qiime tools export \
      --input-path ${params.prefix}_genus_table.qza \
      --output-path genus_export

    biom convert \
      --input-fp genus_export/feature-table.biom \
      --output-fp ${params.prefix}_genus_counts.tsv \
      --to-tsv
    """
}

// --------------------------------------------------
// EXPORT  (merged table + taxonomy to TSV)
// --------------------------------------------------
process export_outputs {

    tag "${params.prefix} - export"

    publishDir params.results_dir, mode: 'copy'

    cpus 1
    memory '4 GB'

    input:
    path table
    path taxonomy

    output:
    path("${params.prefix}_feature-table.tsv")
    path("${params.prefix}_exported-taxonomy.tsv")

    script:
    """
    mkdir -p table_export taxonomy_export

    qiime tools export \
      --input-path $table \
      --output-path table_export

    biom convert \
      --input-fp table_export/feature-table.biom \
      --output-fp ${params.prefix}_feature-table.tsv \
      --to-tsv

    qiime tools export \
      --input-path $taxonomy \
      --output-path taxonomy_export

    cp taxonomy_export/taxonomy.tsv ${params.prefix}_exported-taxonomy.tsv
    """
}

// --------------------------------------------------
// WORKFLOWS
// --------------------------------------------------

// import + demux summary for every batch (for eyeballing quality before denoise)
workflow import_only {

    Channel.fromPath("${params.manifest_dir}/manifest_*.tsv")
        .map { f -> tuple(f.baseName.replaceFirst(/^manifest_/, ''), f) }
        .set { manifest_ch }

    import_reads(manifest_ch)
}


workflow full_run {

    // one tuple (batch_id, manifest) per flowcell manifest
    Channel.fromPath("${params.manifest_dir}/manifest_*.tsv")
        .map { f -> tuple(f.baseName.replaceFirst(/^manifest_/, ''), f) }
        .set { manifest_ch }

    Channel.fromPath(params.classifier).set { classifier_ch }

    // per-batch import + denoise (scatter)
    import_reads(manifest_ch)
    dada2_denoise(import_reads.out.demux_qza)

    // gather all per-batch tables & rep-seqs, then merge (gather)
    merge_batches(
        dada2_denoise.out.table.collect(),
        dada2_denoise.out.repseqs.collect()
    )

    // downstream once on merged artifacts
    summarise_outputs(
        merge_batches.out.table,
        merge_batches.out.repseqs
    )

    assign_taxonomy(
        merge_batches.out.repseqs,
        classifier_ch
    )

    export_outputs(
        merge_batches.out.table,
        assign_taxonomy.out.taxonomy
    )

    collapse_genus(
        merge_batches.out.table,
        assign_taxonomy.out.taxonomy
    )
}
