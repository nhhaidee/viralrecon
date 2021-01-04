process BCFTOOLS_VARIANTS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bcftools", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".txt")) "bcftools_stats/$filename"
                      else filename
                }

    //when:
    //!params.skip_variants && 'bcftools' in callers

    input:
    tuple val(sample), val(single_end), path(bam), path(fasta) /* from ch_markdup_bam_bcftools
    path fasta from ch_fasta */

    output:
    tuple val(sample), val(single_end), path("*.vcf.gz*"), emit: ch_bcftools_variants /* into ch_bcftools_variants_consensus,
                                                               ch_bcftools_variants_snpeff,
                                                               ch_bcftools_variants_intersect */
    path "*.bcftools_stats.txt", emit:  ch_bcftools_variants_mqc

    script:
    """
    echo "$sample" > sample_name.list
    bcftools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth $params.mpileup_depth \\
        --fasta-ref $fasta \\
        --min-BQ $params.min_base_qual \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        ${bam[0]} \\
        | bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${sample}.vcf.gz --output-type z --include 'INFO/DP>=$params.min_coverage'
    tabix -p vcf -f ${sample}.vcf.gz
    bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
    """
}

process BCFTOOLS_CONSENSUS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bcftools/consensus", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "base_qc/$filename"
                      else if (filename.endsWith(".pdf")) "base_qc/$filename"
                      else filename
                }

    //when:
    //!params.skip_variants && 'bcftools' in callers

    input:
    tuple val(sample), val(single_end), path(bam), path(vcf), path(fasta) /*from ch_markdup_bam_bcftools_consensus.join(ch_bcftools_variants_consensus, by: [0,1])
    path fasta from ch_fasta */
    val (index_base)

    output:
    tuple val(sample), val(single_end), path("*consensus.masked.fa"), emit: ch_bcftools_consensus_masked
    path "*.{consensus.fa,tsv,pdf}"

    script:
    """
    cat $fasta | bcftools consensus ${vcf[0]} > ${sample}.consensus.fa
    bedtools genomecov \\
        -bga \\
        -ibam ${bam[0]} \\
        -g $fasta \\
        | awk '\$4 < $params.min_coverage' | bedtools merge > ${sample}.mask.bed
    bedtools maskfasta \\
        -fi ${sample}.consensus.fa \\
        -bed ${sample}.mask.bed \\
        -fo ${sample}.consensus.masked.fa
    sed -i 's/${index_base}/${sample}/g' ${sample}.consensus.masked.fa
    header=\$(head -n1 ${sample}.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sample}/g" ${sample}.consensus.masked.fa
    plot_base_density.r --fasta_files ${sample}.consensus.masked.fa --prefixes $sample --output_dir ./
    """
}

process BCFTOOLS_SNPEFF {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bcftools/snpeff", mode: params.publish_dir_mode

    //when:
    //!params.skip_variants && 'bcftools' in callers && params.gff && !params.skip_snpeff

    input:
    tuple val(sample), val(single_end), path(vcf), file(db), file(config) /*from ch_bcftools_variants_snpeff
    tuple file(db), file(config) from ch_snpeff_db_bcftools*/
    val (index_base)

    output:
    path "*.snpEff.csv", emit: ch_bcftools_snpeff_mqc
    path "*.vcf.gz*"
    path "*.{txt,html}"

    script:
    """
    snpEff ${index_base} \\
        -config $config \\
        -dataDir $db \\
        ${vcf[0]} \\
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
    	"""
}

process BCFTOOLS_QUAST {
    label 'process_medium'
    publishDir "${params.outdir}/variants/bcftools", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (!filename.endsWith(".tsv")) filename
                }

   // when:
    //!params.skip_variants && 'bcftools' in callers && !params.skip_variants_quast

    input:
    path consensus /*from ch_bcftools_consensus_masked.collect{ it[2] } */
    path fasta /*from ch_fasta*/
    path gff /*from ch_gff*/

    output:
    path "quast"
    path "report.tsv", emit: ch_bcftools_quast_mqc

    script:
    features = params.gff ? "--features $gff" : ""
    """
    quast.py \\
        --output-dir quast \\
        -r $fasta \\
        $features \\
        --threads $task.cpus \\
        ${consensus.join(' ')}
    ln -s quast/report.tsv
    """
}

