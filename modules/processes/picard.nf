process PICARD_MARKDUPLICATES {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                        else filename
                }

    when:
    !params.skip_variants

    input:
    tuple val(sample), val(single_end), path(bam), path(fasta) /*from ch_ivar_trim_bam */
    //path fasta /* from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}") , emit: ch_markdup_bam /*into ch_markdup_bam_metrics,
                                                                            ch_markdup_bam_mosdepth_genome,
                                                                            ch_markdup_bam_mosdepth_amplicon,
                                                                            ch_markdup_bam_mpileup,
                                                                            ch_markdup_bam_varscan2_consensus,
                                                                            ch_markdup_bam_bcftools,
                                                                            ch_markdup_bam_bcftools_consensus*/
    path "*.{flagstat,idxstats,stats}", emit: ch_markdup_bam_flagstat_mqc
    path "*.txt", emit: ch_markdup_bam_metrics_mqc

    script:
    def avail_mem = 3
    if (!task.memory) {
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
    } else {
        avail_mem = task.memory.toGiga()
    }
    prefix = params.protocol == 'amplicon' ? "${sample}.trim.mkD" : "${sample}.mkD"
    keep_dup = params.filter_dups ? "true" : "false"
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=${bam[0]} \\
        OUTPUT=${prefix}.sorted.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=$keep_dup \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    samtools index ${prefix}.sorted.bam
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}

process PICARD_METRICS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam/picard_metrics", mode: params.publish_dir_mode

    when:
    !params.skip_variants && !params.skip_picard_metrics

    input:
    tuple val(sample), val(single_end), path(bam), path(fasta) /*from ch_markdup_bam_metrics */
    //path fasta /*from ch_fasta */

    output:
    path "*metrics", emit: ch_picard_metrics_mqc
    path "*.pdf"

    script:
    def avail_mem = 3
    if (!task.memory) {
        log.info "[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
    } else {
        avail_mem = task.memory.toGiga()
    }
    suffix = params.skip_markduplicates ? "" : ".mkD"
    prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}" : "${sample}${suffix}"
    """
    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
        INPUT=${bam[0]} \\
        OUTPUT=${prefix}.CollectMultipleMetrics \\
        REFERENCE_SEQUENCE=$fasta \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    picard -Xmx${avail_mem}g CollectWgsMetrics \\
        COVERAGE_CAP=1000000 \\
        INPUT=${bam[0]} \\
        OUTPUT=${prefix}.CollectWgsMetrics.coverage_metrics \\
        REFERENCE_SEQUENCE=$fasta \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    """
}