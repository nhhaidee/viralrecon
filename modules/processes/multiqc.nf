process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith("assembly_metrics_mqc.tsv")) "assembly/$filename"
                      else if (filename.endsWith("variants_metrics_mqc.tsv")) "variants/$filename"
                      else "multiqc/$filename"
                }

    //when:
    //!params.skip_multiqc

    input:
    val(custom_runName)
    path (multiqc_config) /*from ch_multiqc_config */
    path (mqc_custom_config) /*from ch_multiqc_custom_config.collect().ifEmpty([]) */
    path ('fastqc/*') /*from ch_fastqc_raw_reports_mqc.collect().ifEmpty([]) */
    path ('fastp/log/*') /*from ch_fastp_mqc.collect().ifEmpty([]) */
    path ('fastp/fastqc/*') /*from ch_fastp_fastqc_mqc.collect().ifEmpty([]) */
    path ('bowtie2/log/*') /*from ch_bowtie2_mqc.collect().ifEmpty([]) */
    path ('bowtie2/flagstat/*') /*from ch_sort_bam_flagstat_mqc.collect().ifEmpty([])*/
    path ('ivar/trim/flagstat/*') /*from ch_ivar_trim_flagstat_mqc.collect().ifEmpty([])*/
    path ('ivar/trim/log/*') /*from ch_ivar_trim_log_mqc.collect().ifEmpty([])*/
    path ('picard/markdup/*') /*from ch_markdup_bam_flagstat_mqc.collect().ifEmpty([])*/
    path ('picard/metrics/*') /*from ch_markdup_bam_metrics_mqc.collect().ifEmpty([])*/
    path ('picard/metrics/*') /*from ch_picard_metrics_mqc.collect().ifEmpty([])*/
    path ('mosdepth/genome/*') /*from ch_mosdepth_genome_mqc.collect().ifEmpty([])*/
    path ('varscan2/counts/lowfreq/*') /*from ch_varscan2_log_mqc.collect().ifEmpty([])*/
    path ('varscan2/bcftools/highfreq/*') /*from ch_varscan2_bcftools_highfreq_mqc.collect().ifEmpty([])*/
    path ('varscan2/snpeff/highfreq/*') /*from ch_varscan2_snpeff_highfreq_mqc.collect().ifEmpty([])*/
    path ('varscan2/quast/highfreq/*') /*from ch_varscan2_quast_mqc.collect().ifEmpty([])*/
    path ('ivar/variants/counts/lowfreq/*') /*from ch_ivar_count_mqc.collect().ifEmpty([])*/
    path ('ivar/variants/bcftools/highfreq/*') /*from ch_ivar_bcftools_highfreq_mqc.collect().ifEmpty([])*/
    path ('ivar/variants/snpeff/highfreq/*') /*from ch_ivar_snpeff_highfreq_mqc.collect().ifEmpty([])*/
    path ('ivar/consensus/quast/highfreq/*') /*from ch_ivar_quast_mqc.collect().ifEmpty([])*/
    path ('bcftools/variants/bcftools/*') /*from ch_bcftools_variants_mqc.collect().ifEmpty([])*/
    path ('bcftools/variants/snpeff/*') /*from ch_bcftools_snpeff_mqc.collect().ifEmpty([])*/
    path ('bcftools/consensus/quast/*') /*from ch_bcftools_quast_mqc.collect().ifEmpty([])*/
    path ('cutadapt/log/*')
    path ('cutadapt/fastqc/*') 
    path ('kraken2/*') /*from ch_kraken2_report_mqc.collect().ifEmpty([])*/
    path ('spades/bcftools/*') //from ch_spades_vg_bcftools_mqc.collect().ifEmpty([])
    path ('spades/snpeff/*') //from ch_spades_snpeff_mqc.collect().ifEmpty([])
    path ('spades/quast/*') //from ch_quast_spades_mqc.collect().ifEmpty([])
    path ('metaspades/bcftools/*') //from ch_metaspades_vg_bcftools_mqc.collect().ifEmpty([])
    path ('metaspades/snpeff/*') //from ch_metaspades_snpeff_mqc.collect().ifEmpty([])
    path ('metaspades/quast/*') //from ch_quast_metaspades_mqc.collect().ifEmpty([])
    path ('unicycler/bcftools/*') //from ch_unicycler_vg_bcftools_mqc.collect().ifEmpty([])
    path ('unicycler/snpeff/*') //from ch_unicycler_snpeff_mqc.collect().ifEmpty([])
    path ('unicycler/quast/*') //from ch_quast_unicycler_mqc.collect().ifEmpty([])
    path ('minia/bcftools/*') //from ch_minia_vg_bcftools_mqc.collect().ifEmpty([])
    path ('minia/snpeff/*') //from ch_minia_snpeff_mqc.collect().ifEmpty([])
    path ('minia/quast/*') //from ch_quast_minia_mqc.collect().ifEmpty([])
    path ('software_versions/*') //from ch_software_versions_yaml.collect()
    path workflow_summary //from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    
    output:
    path "*multiqc_report.html",emit: ch_multiqc_report
    path "*_data"
    path "*.tsv"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    multiqc_to_custom_tsv.py
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}