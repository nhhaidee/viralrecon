process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".csv")) filename
                      else null
                }

    output:
    path "software_versions_mqc.yaml", emit: software_versions_yaml
    path "software_versions.csv", emit: software_version_csv

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    parallel-fastq-dump --version > v_parallel_fastq_dump.txt
    fastqc --version > v_fastqc.txt
    fastp --version 2> v_fastp.txt
    bowtie2 --version > v_bowtie2.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    mosdepth --version > v_mosdepth.txt
    picard CollectMultipleMetrics --version &> v_picard.txt || true
    ivar -v > v_ivar.txt
    echo \$(varscan 2>&1) > v_varscan.txt
    bcftools -v > v_bcftools.txt
    snpEff -version > v_snpeff.txt
    echo \$(SnpSift 2>&1) > v_snpsift.txt
    quast.py --version > v_quast.txt
    cutadapt --version > v_cutadapt.txt
    kraken2 --version > v_kraken2.txt
    spades.py --version > v_spades.txt
    unicycler --version > v_unicycler.txt
    minia --version > v_minia.txt
    blastn -version > v_blast.txt
    abacas.pl -v &> v_abacas.txt || true
    plasmidID -v > v_plasmidid.txt  || true
    Bandage --version > v_bandage.txt
    minimap2 --version > v_minimap2.txt
    vg version > v_vg.txt
    echo \$(R --version 2>&1) > v_R.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path output_docs /*from ch_output_docs */
    path images /*from ch_output_docs_images */

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}