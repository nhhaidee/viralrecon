process IVAR_TRIM {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".log")) "log/$filename"
                        else params.skip_markduplicates || params.save_align_intermeds ? filename : null
                }

    when:
    !params.skip_variants

    input:
    tuple val(sample), val(single_end), path(bam) /*from ch_sort_bam*/
    path bed /*from ch_amplicon_bed */

    output:
    tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}"), emit: ch_ivar_trim_bam
    path "*.{flagstat,idxstats,stats}", emit: ch_ivar_trim_flagstat_mqc
    path "*.log", emit: ch_ivar_trim_log_mqc

    script:
    exclude_reads = params.ivar_trim_noprimer ? "" : "-e"
    prefix = "${sample}.trim"
    """
    samtools view -b -F 4 ${bam[0]} > ${sample}.mapped.bam
    samtools index ${sample}.mapped.bam
    ivar trim \\
        -i ${sample}.mapped.bam \\
        -b $bed \\
        -m $params.ivar_trim_min_len \\
        -q $params.ivar_trim_min_qual \\
        -s $params.ivar_trim_window_width \\
        $exclude_reads \\
        -p $prefix > ${prefix}.ivar.log
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $prefix ${prefix}.bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
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
    //!params.skip_variants && 'ivar' in callers

    input:
    tuple val(sample), val(single_end), path(mpileup), path(fasta), path(gff)/* from ch_mpileup_ivar_variants */
    path header /*from ch_ivar_variants_header_mqc */
    /* path fasta from ch_fasta */
    /* path gff from ch_gff */

    output:
    tuple val(sample), val(single_end), path("${prefix}.vcf.gz*"), emit: ch_ivar_highfreq /*into ch_ivar_highfreq_snpeff,
                                                                        ch_ivar_highfreq_intersect */
    tuple val(sample), val(single_end), path("${sample}.vcf.gz*"), emit: ch_ivar_lowfreq /*into ch_ivar_lowfreq_snpeff */
    path "${prefix}.bcftools_stats.txt", emit:  ch_ivar_bcftools_highfreq_mqc
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

process IVAR_CONSENSUS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/ivar/consensus", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "base_qc/$filename"
                      else if (filename.endsWith(".pdf")) "base_qc/$filename"
                      else filename
                }

    //when:
    //!params.skip_variants && 'ivar' in callers

    input:
    tuple val(sample), val(single_end), path(mpileup), path(fasta) /* from ch_mpileup_ivar_consensus
    path fasta from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("*.fa"), emit:  ch_ivar_consensus
    path "*.{txt,tsv,pdf}"

    script:
    prefix = "${sample}.AF${params.max_allele_freq}"
    """
    cat $mpileup | ivar consensus -q $params.min_base_qual -t $params.max_allele_freq -m $params.min_coverage -n N -p ${prefix}.consensus
    header=\$(head -n1 ${prefix}.consensus.fa | sed 's/>//g')
    sed -i "s/\${header}/${sample}/g" ${prefix}.consensus.fa
    plot_base_density.r --fasta_files ${prefix}.consensus.fa --prefixes $prefix --output_dir ./
    """
}

process IVAR_SNPEFF {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/ivar/snpeff", mode: params.publish_dir_mode

    //when:
    //!params.skip_variants && 'ivar' in callers && params.gff && !params.skip_snpeff

    input:
    tuple val(sample), val(single_end), path(highfreq_vcf), path(lowfreq_vcf),file(db), file(config)
    val (index_base) /*from ch_ivar_highfreq_snpeff.join(ch_ivar_lowfreq_snpeff, by: [0,1])
    
    tuple file(db), file(config) from ch_snpeff_db_ivar */

    output:
    path "${prefix}.snpEff.csv", emit: ch_ivar_snpeff_highfreq_mqc
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
process IVAR_QUAST {
    label 'process_medium'
    publishDir "${params.outdir}/variants/ivar/quast", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (!filename.endsWith(".tsv")) filename
                }

    //when:
    //!params.skip_variants && 'ivar' in callers && !params.skip_variants_quast

    input:
    path consensus /*from ch_ivar_consensus.collect{ it[2] } */
    path fasta /*from ch_fasta */
    path gff /*from ch_gff */

    output:
    path "AF${params.max_allele_freq}"
    path "report.tsv", emit: ch_ivar_quast_mqc

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