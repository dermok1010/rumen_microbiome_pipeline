// main.nf
nextflow.enable.dsl=2

// --------------------------------------------------
// IMPORT READS
// --------------------------------------------------
process import_reads {

    tag "${params.prefix} - import"

    publishDir params.output_dir, mode: 'copy'

    cpus 2
    memory '8 GB'
    time '1h'

    input:
    path manifest_file

    output:
    path("${params.prefix}_paired_end_demux.qza"), emit: demux_qza
    path("${params.prefix}_paired_end_demux.qzv")

    script:
    """
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path $manifest_file \
      --output-path ${params.prefix}_paired_end_demux.qza \
      --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
      --i-data ${params.prefix}_paired_end_demux.qza \
      --o-visualization ${params.prefix}_paired_end_demux.qzv
    """
}

// --------------------------------------------------
// DADA2
// --------------------------------------------------
process dada2_denoise {

    tag "${params.prefix} - dada2"

    publishDir params.output_dir, mode: 'copy'

    cpus 8
    memory '24 GB'
    time '8h'

    input:
    path demux_qza

    output:
    path("${params.prefix}_table.qza"), emit: table
    path("${params.prefix}_rep_seqs.qza"), emit: repseqs
    path("${params.prefix}_denoising_stats.qza"), emit: stats

    script:
    """
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $demux_qza \
      --p-trim-left-f ${params.trim_left_f} \
      --p-trim-left-r ${params.trim_left_r} \
      --p-trunc-len-f ${params.trunc_len_f} \
      --p-trunc-len-r ${params.trunc_len_r} \
      --p-n-threads ${task.cpus} \
      --o-table ${params.prefix}_table.qza \
      --o-representative-sequences ${params.prefix}_rep_seqs.qza \
      --o-denoising-stats ${params.prefix}_denoising_stats.qza \
      --verbose
    """
}

// --------------------------------------------------
// SUMMARIES
// --------------------------------------------------
process summarise_outputs {

    tag "${params.prefix} - summary"

    publishDir params.output_dir, mode: 'copy'

    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path table
    path repseqs
    path stats_file

    output:
    path("${params.prefix}_ct_table.qzv"), emit: summary_table
    path("${params.prefix}_ct_rep_seqs.qzv"), emit: summary_repseqs
    path("${params.prefix}_ct_stats.qzv"), emit: summary_stats

    script:
    """
    qiime feature-table summarize \
      --i-table $table \
      --o-visualization ${params.prefix}_ct_table.qzv

    qiime feature-table tabulate-seqs \
      --i-data $repseqs \
      --o-visualization ${params.prefix}_ct_rep_seqs.qzv

    qiime metadata tabulate \
      --m-input-file $stats_file \
      --o-visualization ${params.prefix}_ct_stats.qzv
    """
}

// --------------------------------------------------
// TAXONOMY
// --------------------------------------------------
process assign_taxonomy {

    tag "${params.prefix} - taxonomy"

    publishDir params.output_dir, mode: 'copy'

    cpus 4
    memory '16 GB'
    time '4h'

    input:
    path repseqs
    path classifier_file

    output:
    path("${params.prefix}_ct_taxonomy.qza"), emit: taxonomy
    path("${params.prefix}_ct_taxonomy.qzv")

    script:
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier $classifier_file \
      --i-reads $repseqs \
      --o-classification ${params.prefix}_ct_taxonomy.qza \
      --p-n-jobs ${task.cpus}

    qiime metadata tabulate \
      --m-input-file ${params.prefix}_ct_taxonomy.qza \
      --o-visualization ${params.prefix}_ct_taxonomy.qzv
    """
}

// --------------------------------------------------
// EXPORT
// --------------------------------------------------
process export_outputs {

    tag "${params.prefix} - export"

    publishDir params.results_dir, mode: 'copy'

    cpus 1
    memory '4 GB'
    time '1h'

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

workflow import_only {

    Channel.fromPath(params.manifest).set { manifest_ch }

    import_reads(manifest_ch)
}


workflow full_run {

    Channel.fromPath(params.manifest).set { manifest_ch }
    Channel.fromPath(params.classifier).set { classifier_ch }

    import_reads(manifest_ch)

    dada2_denoise(import_reads.out.demux_qza)

    summarise_outputs(
        dada2_denoise.out.table,
        dada2_denoise.out.repseqs,
        dada2_denoise.out.stats
    )

    assign_taxonomy(
        dada2_denoise.out.repseqs,
        classifier_ch
    )

    export_outputs(
        dada2_denoise.out.table,
        assign_taxonomy.out.taxonomy
    )
}
