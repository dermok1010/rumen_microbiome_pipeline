// main.nf
nextflow.enable.dsl=2

// PROCESSES
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
    source activate qiime2-amplicon-2024.2

    qiime tools import \\
      --type 'SampleData[PairedEndSequencesWithQuality]' \\
      --input-path $manifest_file \\
      --output-path ${params.prefix}_paired_end_demux.qza \\
      --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \\
      --i-data ${params.prefix}_paired_end_demux.qza \\
      --o-visualization ${params.prefix}_paired_end_demux.qzv
    """
}

process dada2_denoise {
    tag "${params.prefix} - dada2"
    publishDir params.output_dir, mode: 'copy'
    cpus 8
    memory '8 GB'
    time '2h'

    input:
    path demux_qza

    output:
    path("${params.prefix}_table.qza"), emit: table
    path("${params.prefix}_rep_seqs.qza"), emit: repseqs
    path("${params.prefix}_denoising_stats.qza"), emit: stats

    script:
    """
    source activate qiime2-amplicon-2024.2

    qiime dada2 denoise-paired \\
      --i-demultiplexed-seqs $demux_qza \\
      --p-trim-left-f 17 \\
      --p-trim-left-r 21 \\
      --p-trunc-len-f 230 \\
      --p-trunc-len-r 200 \\
      --p-n-threads 8 \\
      --o-table ${params.prefix}_table.qza \\
      --o-representative-sequences ${params.prefix}_rep_seqs.qza \\
      --o-denoising-stats ${params.prefix}_denoising_stats.qza \\
      --verbose
    """
}

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
    source activate qiime2-amplicon-2024.2

    qiime feature-table summarize \\
      --i-table $table \\
      --o-visualization ${params.prefix}_ct_table.qzv

    qiime feature-table tabulate-seqs \\
      --i-data $repseqs \\
      --o-visualization ${params.prefix}_ct_rep_seqs.qzv

    qiime metadata tabulate \\
      --m-input-file $stats_file \\
      --o-visualization ${params.prefix}_ct_stats.qzv
    """
}

process assign_taxonomy {
    tag "${params.prefix} - taxonomy"
    publishDir params.output_dir, mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    path repseqs
    path classifier_file

    output:
    path("${params.prefix}_ct_taxonomy.qza"), emit: taxonomy
    path("${params.prefix}_ct_taxonomy.qzv")

    script:
    """
    source activate qiime2-amplicon-2024.2

    qiime feature-classifier classify-sklearn \\
      --i-classifier $classifier_file \\
      --i-reads $repseqs \\
      --o-classification ${params.prefix}_ct_taxonomy.qza \\
      --p-n-jobs 4

    qiime metadata tabulate \\
      --m-input-file ${params.prefix}_ct_taxonomy.qza \\
      --o-visualization ${params.prefix}_ct_taxonomy.qzv
    """
}

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
    source activate qiime2-amplicon-2024.2

    qiime tools export \\
      --input-path $table \\
      --output-path .

    biom convert \\
      --input-fp feature-table.biom \\
      --output-fp ${params.prefix}_feature-table.tsv \\
      --to-tsv

    qiime tools export \\
      --input-path $taxonomy \\
      --output-path .

    mv taxonomy.tsv ${params.prefix}_exported-taxonomy.tsv
    """
}

workflow {
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

