
process KRAKEN2_BUILD {
    tag "$db"
    label 'process_high'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    //when:
    //!params.skip_assembly

    output:
    path "$db", emit: ch_kraken2_db

    script:
    db = "kraken2_${params.kraken2_db_name}"
    ftp = params.kraken2_use_ftp ? "--use-ftp" : ""
    """
    kraken2-build --db $db --threads $task.cpus $ftp --download-taxonomy
    kraken2-build --db $db --threads $task.cpus $ftp --download-library $params.kraken2_db_name
    kraken2-build --db $db --threads $task.cpus $ftp --build
    """
}

process KRAKEN2 {
    tag "$db"
    label 'process_high'
    publishDir "${params.outdir}/assembly/kraken2", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if (filename.endsWith(".txt")) filename
                        else params.save_kraken2_fastq ? filename : null
                }

    input:
    tuple val(sample), val(single_end), path(reads), path(db) /* from ch_fastp_kraken2 */
    /*path db from ch_kraken2_db */

    output:
    tuple val(sample), val(single_end), path("*.viral*"), emit: ch_kraken2 /* into ch_kraken2_spades,
                                                                ch_kraken2_metaspades,
                                                                ch_kraken2_unicycler,
                                                                ch_kraken2_minia */
    path "*.report.txt", emit: ch_kraken2_report_mqc
    path "*.host*"


    script:
    pe = single_end ? "" : "--paired"
    classified = single_end ? "${sample}.host.fastq" : "${sample}.host#.fastq"
    unclassified = single_end ? "${sample}.viral.fastq" : "${sample}.viral#.fastq"
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --unclassified-out $unclassified \\
        --classified-out $classified \\
        --report ${sample}.kraken2.report.txt \\
        --report-zero-counts \\
        $pe \\
        --gzip-compressed \\
        $reads
    pigz -p $task.cpus *.fastq
    """
}