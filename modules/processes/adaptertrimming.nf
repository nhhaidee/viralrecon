process FASTP {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/preprocess/fastp", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".json")) filename
                    else if (filename.endsWith(".fastp.html")) filename
                    else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
                    else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                    else if (filename.endsWith(".log")) "log/$filename"
                    else params.save_trimmed ? filename : null
                }

    when:
    !params.skip_variants || !params.skip_assembly

    input:
    tuple val(sample), val(single_end), path(reads)/* from ch_cat_fastp */

    output:
    tuple val(sample), val(single_end), path("*.trim.fastq.gz"), emit: ch_fastp /* into ch_fastp_bowtie2,
                                                                    ch_fastp_cutadapt,
                                                                    ch_fastp_kraken2 */
    path "*.json", emit: ch_fastp_mqc
    path "*_fastqc.{zip,html}", emit: ch_fastp_fastqc_mqc
    path "*.{log,fastp.html}"
    path "*.fail.fastq.gz"

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    autodetect = single_end ? "" : "--detect_adapter_for_pe"
    """
    IN_READS='--in1 ${sample}.fastq.gz'
    OUT_READS='--out1 ${sample}.trim.fastq.gz --failed_out ${sample}.fail.fastq.gz'
    if $single_end; then
        [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
    else
        [ ! -f  ${sample}_1.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.fastq.gz
        [ ! -f  ${sample}_2.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.fastq.gz
        IN_READS='--in1 ${sample}_1.fastq.gz --in2 ${sample}_2.fastq.gz'
        OUT_READS='--out1 ${sample}_1.trim.fastq.gz --out2 ${sample}_2.trim.fastq.gz --unpaired1 ${sample}_1.fail.fastq.gz --unpaired2 ${sample}_2.fail.fastq.gz'
    fi
    fastp \\
        \$IN_READS \\
        \$OUT_READS \\
        $autodetect \\
        --cut_front \\
        --cut_tail \\
        --cut_mean_quality $params.cut_mean_quality \\
        --qualified_quality_phred $params.qualified_quality_phred \\
        --unqualified_percent_limit $params.unqualified_percent_limit \\
        --length_required $params.min_trim_length \\
        --trim_poly_x \\
        --thread $task.cpus \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html \\
        2> ${sample}.fastp.log
    fastqc --quiet --threads $task.cpus *.trim.fastq.gz
    """
}