process MOSDEPTH_GENOME {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam/mosdepth/genome", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".pdf")) "plots/$filename"
                      else if (filename.endsWith(".tsv")) "plots/$filename"
                      else filename
                }

    when:
    !params.skip_variants && !params.skip_mosdepth

    input:
    tuple val(sample), val(single_end), path(bam) /*from ch_markdup_bam_mosdepth_genome */

    output:
    path "*.global.dist.txt", emit: ch_mosdepth_genome_mqc
    path "*.{txt,gz,csi,tsv,pdf}"

    script:
    suffix = params.skip_markduplicates ? "" : ".mkD"
    prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}.genome" : "${sample}${suffix}.genome"
    plot_suffix = params.protocol == 'amplicon' ? ".trim${suffix}.genome" : "${suffix}.genome"
    """
    mosdepth \\
        --by 200 \\
        --fast-mode \\
        $prefix \\
        ${bam[0]}
    plot_mosdepth_regions.r \\
        --input_files ${prefix}.regions.bed.gz \\
        --input_suffix ${plot_suffix}.regions.bed.gz \\
        --output_dir ./ \\
        --output_suffix ${plot_suffix}.regions
    """
}

process MOSDEPTH_AMPLICON {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam/mosdepth/amplicon", mode: params.publish_dir_mode

    when:
    !params.skip_variants && !params.skip_mosdepth

    input:
    tuple val(sample), val(single_end), path(bam) /* from ch_markdup_bam_mosdepth_amplicon */
    path bed/* from ch_amplicon_bed */

    output:
    path "*.regions.bed.gz", emit: ch_mosdepth_amplicon_region_bed
    path "*.{txt,gz,csi}"

    script:
    suffix = params.skip_markduplicates ? "" : ".mkD"
    prefix = "${sample}.trim${suffix}.amplicon"
    """
    collapse_amplicon_bed.py \\
        --left_primer_suffix $params.amplicon_left_suffix \\
        --right_primer_suffix $params.amplicon_right_suffix \\
        $bed \\
        amplicon.collapsed.bed
    mosdepth \\
        --by amplicon.collapsed.bed \\
        --fast-mode \\
        --use-median \\
        --thresholds 0,1,10,50,100,500 \\
        ${prefix} \\
        ${bam[0]}
    """
}

process MOSDEPTH_AMPLICON_PLOT {
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam/mosdepth/amplicon/plots", mode: params.publish_dir_mode

    when:
    !params.skip_variants && !params.skip_mosdepth

    input:
    path bed /* from ch_mosdepth_amplicon_region_bed.collect() */

    output:
    path "*.{tsv,pdf}"

    script:
    suffix = params.skip_markduplicates ? "" : ".mkD"
    suffix = ".trim${suffix}.amplicon"
    """
    plot_mosdepth_regions.r \\
        --input_files ${bed.join(',')} \\
        --input_suffix ${suffix}.regions.bed.gz \\
        --output_dir ./ \\
        --output_suffix ${suffix}.regions
    """
}
