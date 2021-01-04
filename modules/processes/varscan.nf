process VARSCAN2 {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/varscan2", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".log")) "log/$filename"
                      else if (filename.endsWith(".txt")) "bcftools_stats/$filename"
                      else filename
                }

    //when:
    //!params.skip_variants && 'varscan2' in callers

    input:
    tuple val(sample), val(single_end), path(mpileup), path(fasta) /*from ch_mpileup_varscan2 */
    /*path fasta from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("${prefix}.vcf.gz*"), emit: ch_varscan2_highfreq /* into ch_varscan2_highfreq_consensus,
                                                                       ch_varscan2_highfreq_snpeff,
                                                                       ch_varscan2_highfreq_intersect */
    tuple val(sample), val(single_end), path("${sample}.vcf.gz*"), emit: ch_varscan2_lowfreq
    path "${prefix}.bcftools_stats.txt", emit: ch_varscan2_bcftools_highfreq_mqc
    path "*.varscan2.log", emit: ch_varscan2_log_mqc
    path "${sample}.bcftools_stats.txt"

    script:
    prefix = "${sample}.AF${params.max_allele_freq}"
    strand = params.protocol != 'amplicon' && params.varscan2_strand_filter ? "--strand-filter 1" : "--strand-filter 0"
    """
    echo "$sample" > sample_name.list
    varscan mpileup2cns \\
        $mpileup \\
        --min-coverage $params.min_coverage \\
        --min-reads2 5 \\
        --min-avg-qual $params.min_base_qual \\
        --min-var-freq $params.min_allele_freq \\
        --p-value 0.99 \\
        --output-vcf 1 \\
        --vcf-sample-list sample_name.list \\
        --variants \\
        $strand \\
        2> ${sample}.varscan2.log \\
        | bgzip -c > ${sample}.vcf.gz
    tabix -p vcf -f ${sample}.vcf.gz
    bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
    sed -i.bak '/LC_ALL/d' ${sample}.varscan2.log
    bcftools filter \\
        -i 'FORMAT/AD / (FORMAT/AD + FORMAT/RD) >= $params.max_allele_freq' \\
        --output-type z \\
        --output ${prefix}.vcf.gz \\
        ${sample}.vcf.gz
    tabix -p vcf -f ${prefix}.vcf.gz
    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
    """
}

process VARSCAN2_CONSENSUS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/varscan2/consensus", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "base_qc/$filename"
                      else if (filename.endsWith(".pdf")) "base_qc/$filename"
                      else filename
                }

    //when:
    //!params.skip_variants && 'varscan2' in callers

    input:
    tuple val(sample), val(single_end), path(bam), path(vcf), path(fasta) /*from ch_markdup_bam_varscan2_consensus.join(ch_varscan2_highfreq_consensus, by: [0,1])
    path fasta from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("*consensus.masked.fa"), emit: ch_varscan2_consensus
    path "*.{consensus.fa,tsv,pdf}"

    script:
    prefix = "${sample}.AF${params.max_allele_freq}"
    """
    cat $fasta | bcftools consensus ${vcf[0]} > ${prefix}.consensus.fa
    bedtools genomecov \\
        -bga \\
        -ibam ${bam[0]} \\
        -g $fasta \\
        | awk '\$4 < $params.min_coverage' | bedtools merge > ${prefix}.mask.bed
    bedtools maskfasta \\
        -fi ${prefix}.consensus.fa \\
        -bed ${prefix}.mask.bed \\
        -fo ${prefix}.consensus.masked.fa
    header=\$(head -n 1 ${prefix}.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sample}/g" ${prefix}.consensus.masked.fa
    plot_base_density.r --fasta_files ${prefix}.consensus.masked.fa --prefixes $prefix --output_dir ./
    """
}

process VARSCAN2_SNPEFF {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/varscan2/snpeff", mode: params.publish_dir_mode

    //when:
    //!params.skip_variants && 'varscan2' in callers && params.gff && !params.skip_snpeff

    input:
    tuple val(sample), val(single_end), path(highfreq_vcf), path(lowfreq_vcf), file(db), file(config)/* from ch_varscan2_highfreq_snpeff.join(ch_varscan2_lowfreq_snpeff, by: [0,1]) */
     /*tuple file(db), file(config) from ch_snpeff_db_varscan2 */
    val (index_base)

    output:
    path "${prefix}.snpEff.csv", emit: ch_varscan2_snpeff_highfreq_mqc
    path "${sample}.snpEff.csv"
    path "*.vcf.gz*"
    path "*.{txt,html}"

    script:
    prefix = "${sample}.AF${params.max_allele_freq}"
    """
    snpEff ${index_base} \\
        -config $config \\
        -dataDir $db \\
        ${lowfreq_vcf[0]} \\
        -csvStats ${sample}.snpEff.csv \\
        | bgzip -c > ${sample}.snpEff.vcf.gz
    tabix -p vcf -f ${sample}.snpEff.vcf.gz
    mv snpEff_summary.html ${sample}.snpEff.summary.html
    SnpSift extractFields -s "," \\
        -e "." \\
        ${sample}.snpEff.vcf.gz \\
        CHROM POS REF ALT \\
        "ANN[*].GENE" "ANN[*].GENEID" \\
        "ANN[*].IMPACT" "ANN[*].EFFECT" \\
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
        > ${sample}.snpSift.table.txt
    snpEff ${index_base} \\
        -config $config \\
        -dataDir $db \\
        ${highfreq_vcf[0]} \\
        -csvStats ${prefix}.snpEff.csv \\
        | bgzip -c > ${prefix}.snpEff.vcf.gz
    tabix -p vcf -f ${prefix}.snpEff.vcf.gz
    mv snpEff_summary.html ${prefix}.snpEff.summary.html
    SnpSift extractFields -s "," \\
        -e "." \\
        ${prefix}.snpEff.vcf.gz \\
        CHROM POS REF ALT \\
        "ANN[*].GENE" "ANN[*].GENEID" \\
        "ANN[*].IMPACT" "ANN[*].EFFECT" \\
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
        > ${prefix}.snpSift.table.txt
    	"""
}

process VARSCAN2_QUAST {
    label 'process_medium'
    publishDir "${params.outdir}/variants/varscan2/quast", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (!filename.endsWith(".tsv")) filename
                }

    //when:
   // !params.skip_variants && 'varscan2' in callers && !params.skip_variants_quast

    input:
    path consensus/* from ch_varscan2_consensus.collect{ it[2] } */
    path fasta /* from ch_fasta */
    path gff/* from ch_gff */

    output:
    path "AF${params.max_allele_freq}"
    path "report.tsv", emit: ch_varscan2_quast_mqc

    script:
    features = params.gff ? "--features $gff" : ""
    """
    quast.py \\
        --output-dir AF${params.max_allele_freq} \\
        -r $fasta \\
        $features \\
        --threads $task.cpus \\
        ${consensus.join(' ')}
    ln -s AF${params.max_allele_freq}/report.tsv
    """
}

process IVAR_VARIANTS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/ivar", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".bcftools_stats.txt")) "bcftools_stats/$filename"
                      else if (filename.endsWith(".log")) "log/$filename"
                      else if (filename.endsWith("_mqc.tsv")) null
                      else filename
                }

    //when:
   // !params.skip_variants && 'ivar' in callers

    input:
    tuple val(sample), val(single_end), path(mpileup) /*from ch_mpileup_ivar_variants */
    path header /*from ch_ivar_variants_header_mqc */
    path fasta /*from ch_fasta */
    path gff /*from ch_gff */

    output:
    tuple val(sample), val(single_end), path("${prefix}.vcf.gz*"), emit: ch_ivar_highfreq /*into ch_ivar_highfreq_snpeff,
                                                                       ch_ivar_highfreq_intersect*/
    tuple val(sample), val(single_end), path("${sample}.vcf.gz*"), emit: ch_ivar_lowfreq /* into ch_ivar_lowfreq_snpeff*/
    path "${prefix}.bcftools_stats.txt", emit: ch_ivar_bcftools_highfreq_mqc
    path "${sample}.variant.counts_mqc.tsv", emit: ch_ivar_count_mqc
    path "${sample}.bcftools_stats.txt"
    path "${sample}.tsv"
    path "*.log"

    script:
    features = params.gff ? "-g $gff" : ""
    prefix = "${sample}.AF${params.max_allele_freq}"
    """
    cat $mpileup | ivar variants -q $params.min_base_qual -t $params.min_allele_freq -m $params.min_coverage -r $fasta $features -p $sample
    ivar_variants_to_vcf.py ${sample}.tsv ${sample}.vcf > ${sample}.variant.counts.log
    bgzip -c ${sample}.vcf > ${sample}.vcf.gz
    tabix -p vcf -f ${sample}.vcf.gz
    bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
    cat $header ${sample}.variant.counts.log > ${sample}.variant.counts_mqc.tsv
    ivar_variants_to_vcf.py ${sample}.tsv ${prefix}.vcf --pass_only --allele_freq_thresh $params.max_allele_freq > ${prefix}.variant.counts.log
    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf -f ${prefix}.vcf.gz
    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
    """
}