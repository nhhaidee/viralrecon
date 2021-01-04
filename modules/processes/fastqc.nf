process CAT_FASTQ {
    tag "$sample"

    input:
    tuple val(sample), val(single_end), path(reads)

    output:
    tuple val(sample), val(single_end), path("*.merged.fastq.gz"), emit: ch_cat_fastq

    script:
    readList = reads.collect{it.toString()}
    if (!single_end) {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${sample}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${sample}_2.merged.fastq.gz
            """
        } else {
            """
            ln -s ${reads[0]} ${sample}_1.merged.fastq.gz
            ln -s ${reads[1]} ${sample}_2.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 1) {
            """
            cat ${readList.sort().join(' ')} > ${sample}.merged.fastq.gz
            """
        } else {
            """
            ln -s $reads ${sample}.merged.fastq.gz
            """
        }
    }
}

process FASTQC {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/preprocess/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.endsWith(".zip") ? "zips/$filename" : filename
                }

    when:
    !params.skip_fastqc

    input:
    tuple val(sample), val(single_end), path(reads)

    output:
    path "*.{zip,html}", emit: ch_fastqc_raw_reports_mqc

    script:
    """
    fastqc --quiet --threads $task.cpus *.fastq.gz
    """
}