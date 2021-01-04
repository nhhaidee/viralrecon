process SRA_FASTQ_FTP {
    tag "$sample"
    label 'process_medium'
    label 'error_retry'
    publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if (filename.endsWith(".md5")) "md5/$filename"
                        else params.save_sra_fastq ? filename : null
                }

    when:
    is_ftp

    input:
    tuple val(sample), val(single_end), val(is_sra), val(is_ftp), val(fastq), val(md5) /*from ch_reads_sra_ftp*/

    output:
    tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") /*into ch_sra_fastq_ftp*/
    path "*.md5"

    script:
    if (single_end) {
        """
        curl -L ${fastq[0]} -o ${sample}.fastq.gz
        echo "${md5[0]}  ${sample}.fastq.gz" > ${sample}.fastq.gz.md5
        md5sum -c ${sample}.fastq.gz.md5
        """
    } else {
        """
        curl -L ${fastq[0]} -o ${sample}_1.fastq.gz
        echo "${md5[0]}  ${sample}_1.fastq.gz" > ${sample}_1.fastq.gz.md5
        md5sum -c ${sample}_1.fastq.gz.md5
        curl -L ${fastq[1]} -o ${sample}_2.fastq.gz
        echo "${md5[1]}  ${sample}_2.fastq.gz" > ${sample}_2.fastq.gz.md5
        md5sum -c ${sample}_2.fastq.gz.md5
        """
    }
}

process SRA_FASTQ_DUMP {
    tag "$sample"
    label 'process_medium'
    label 'error_retry'
    publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if (filename.endsWith(".log")) "log/$filename"
                        else params.save_sra_fastq ? filename : null
                }

    when:
    !is_ftp

    input:
    tuple val(sample), val(single_end), val(is_sra), val(is_ftp) /*from ch_reads_sra_dump.map { it[0..3] }*/

    output:
    tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") /*into ch_sra_fastq_dump*/
    path "*.log"

    script:
    prefix = "${sample.split('_')[0..-2].join('_')}"
    pe = single_end ? "" : "--readids --split-e"
    rm_orphan = single_end ? "" : "[ -f  ${prefix}.fastq.gz ] && rm ${prefix}.fastq.gz"
    """
    parallel-fastq-dump \\
        --sra-id $prefix \\
        --threads $task.cpus \\
        --outdir ./ \\
        --tmpdir ./ \\
        --gzip \\
        $pe \\
        > ${prefix}.fastq_dump.log
    $rm_orphan
    """
}

