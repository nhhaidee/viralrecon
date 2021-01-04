process SORT_BAM {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                      else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                      else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                      else (params.protocol != 'amplicon' && params.skip_markduplicates) || params.save_align_intermeds ? filename : null
                }

    when:
    !params.skip_variants

    input:
    tuple val(sample), val(single_end), path(bam) /*from ch_bowtie2_bam*/

    output:
    tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}"), path("*.flagstat"), emit: ch_sort_bam
    path "*.{flagstat,idxstats,stats}", emit: ch_sort_bam_flagstat_mqc

    script:
    """
    samtools sort -@ $task.cpus -o ${sample}.sorted.bam -T $sample $bam
    samtools index ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}.sorted.bam.flagstat
    samtools idxstats ${sample}.sorted.bam > ${sample}.sorted.bam.idxstats
    samtools stats ${sample}.sorted.bam > ${sample}.sorted.bam.stats
    """
}

process SAMTOOLS_MPILEUP {
    tag "$sample"
    label 'process_medium'
    if (params.save_mpileup) {
        publishDir "${params.outdir}/variants/bam/mpileup", mode: params.publish_dir_mode
    }

    when:
    !params.skip_variants

    input:
    tuple val(sample), val(single_end), path(bam), path(fasta) /*from ch_markdup_bam_mpileup */
    /*path fasta from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("*.mpileup"), emit:  ch_mpileup /*into ch_mpileup_varscan2,
                                                               ch_mpileup_ivar_variants,
                                                               ch_mpileup_ivar_consensus,
                                                               ch_mpileup_ivar_bcftools */

    script:
    suffix = params.skip_markduplicates ? "" : ".mkD"
    prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}" : "${sample}${suffix}"
    """
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth $params.mpileup_depth \\
        --fasta-ref $fasta \\
        --min-BQ $params.min_base_qual \\
        --output ${prefix}.mpileup \\
        ${bam[0]}
    """
}