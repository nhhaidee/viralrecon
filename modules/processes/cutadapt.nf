process CUTADAPT {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/assembly/cutadapt", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".html")) "fastqc/$filename"
                    else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                    else if (filename.endsWith(".log")) "log/$filename"
                    else params.save_trimmed ? filename : null
                }

    input:
    tuple val(sample), val(single_end), path(reads)/* from ch_fastp_cutadapt */
    path amplicons/* from ch_amplicon_fasta */

    output:
    tuple val(sample), val(single_end), path("*.ptrim.fastq.gz"), emit: ch_cutadapt_kraken2
    path "*.{zip,html}", emit: ch_cutadapt_fastqc_mqc
    path "*.log", emit: ch_cutadapt_mqc

    script:
    adapters = single_end ? "-a file:primers.fasta" : "-a file:primers.fasta -A file:primers.fasta"
    out_reads = single_end ? "-o ${sample}.ptrim.fastq.gz" : "-o ${sample}_1.ptrim.fastq.gz -p ${sample}_2.ptrim.fastq.gz"
    """
    sed -r '/^[ACTGactg]+\$/ s/\$/X/g' $amplicons > primers.fasta
    cutadapt \\
        --cores $task.cpus \\
        --overlap 5 \\
        --minimum-length 30 \\
        --error-rate 0.1 \\
        $adapters \\
        $out_reads \\
        $reads \\
        > ${sample}.cutadapt.log
    fastqc --quiet --threads $task.cpus *.ptrim.fastq.gz
    """
}