process BCFTOOLS_ISEC {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/variants/intersect", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(single_end), path('varscan2/*'), path('ivar/*'), path('bcftools/*') /*from ch_varscan2_highfreq_intersect*/

    output:
    path "$sample"

    script:
    """
    bcftools isec  \\
        --nfiles +2 \\
        --output-type z \\
        -p $sample \\
        */*.vcf.gz
    """
}