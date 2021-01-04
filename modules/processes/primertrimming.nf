process FGBIO_TRIM_PRIMERS {
    
    echo true
    tag "$sample"

    input:
    path primer_info_tab/* from ch_primer_info*/
    tuple val(sample), val(single_end), path(bam) /*from ch_bowtie2_bam*/

    output:
    tuple val(sample), val(single_end), path("*.primerTrim.bam"), emit: ch_fgbio_trim_primers_bam

    script:
    """
    fgbio --sam-validation-stringency=LENIENT TrimPrimers -i $bam -o ${sample}.primerTrim.bam -p $primer_info_tab -H true
    """
}