process METASPADES {
    tag "$sample"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades", mode: params.publish_dir_mode,
    saveAs: { filename ->
                  if (filename.endsWith(".png")) "bandage/$filename"
                  else if (filename.endsWith(".svg")) "bandage/$filename"
                  else filename
            }

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end

    input:
    tuple val(sample), val(single_end), path(reads) //from ch_kraken2_metaspades

    output:
    tuple val(sample), val(single_end), path("*scaffolds.fa"), emit: ch_metaspades /* into ch_metaspades_blast,
                                                                   ch_metaspades_abacas,
                                                                   ch_metaspades_plasmidid,
                                                                   ch_metaspades_quast,
                                                                   ch_metaspades_vg */
    path "*assembly.{gfa,png,svg}"


    script:
    """
    spades.py \\
        --meta \\
        --threads $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ./
    mv scaffolds.fasta ${sample}.scaffolds.fa
    mv assembly_graph_with_scaffolds.gfa ${sample}.assembly.gfa
    if [ -s ${sample}.assembly.gfa ]
    then
        Bandage image ${sample}.assembly.gfa ${sample}.assembly.png --height 1000
        Bandage image ${sample}.assembly.gfa ${sample}.assembly.svg --height 1000
    fi
    """
}
process METASPADES_BLAST {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades/blast", mode: params.publish_dir_mode

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_blast

    input:
    tuple val(sample), val(single_end), path(scaffold), path (db)// from ch_metaspades_blast
    //path db from ch_blast_db
    path header //from ch_blast_outfmt6_header
    val (fasta_base)

    output:
    path "*.blast*"

    script:
    """
    blastn \\
        -num_threads $task.cpus \\
        -db $db/$fasta_base \\
        -query $scaffold \\
        -outfmt \'6 stitle std slen qlen qcovs\' \\
        -out ${sample}.blast.txt
    awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' ${sample}.blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > ${sample}.blast.filt.txt
    cat $header ${sample}.blast.filt.txt > ${sample}.blast.filt.header.txt
    """
}

process METASPADES_ABACAS {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades/abacas", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
                      else filename
                }

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_abacas

    input:
    tuple val(sample), val(single_end), path(scaffold), path (fasta)  //from ch_metaspades_abacas
    //path fasta from ch_fasta

    output:
    path "*.abacas*"

    script:
    """
    abacas.pl -r $fasta -q $scaffold -m -p nucmer -o ${sample}.abacas
    mv nucmer.delta ${sample}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${sample}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${sample}.abacas.nucmer.tiling
    mv unused_contigs.out ${sample}.abacas.unused.contigs.out
    """
}

process METASPADES_PLASMIDID {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades/plasmidid", mode: params.publish_dir_mode

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_plasmidid

    input:
    tuple val(sample), val(single_end), path(scaffold), path (fasta) //from ch_metaspades_plasmidid.filter { it.size() > 0 }
   // path fasta from ch_fasta

    output:
    path "$sample"

    script:
    """
    plasmidID -d $fasta -s $sample -c $scaffold --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o .
    mv NO_GROUP/$sample ./$sample
    """
}

process METASPADES_QUAST {
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (!filename.endsWith(".tsv")) filename
                }

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_assembly_quast

    input:
    path scaffolds //from ch_metaspades_quast.collect{ it[2] }
    path fasta //from ch_fasta
    path gff //from ch_gff

    output:
    path "quast"
    path "report.tsv", emit:  ch_quast_metaspades_mqc

    script:
    features = params.gff ? "--features $gff" : ""
    """
    quast.py \\
        --output-dir quast \\
        -r $fasta \\
        $features \\
        --threads $task.cpus \\
        ${scaffolds.join(' ')}
    ln -s quast/report.tsv
    """
}

process METASPADES_VG {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades/variants", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".txt")) "bcftools_stats/$filename"
                      else if (filename.endsWith(".png")) "bandage/$filename"
                      else if (filename.endsWith(".svg")) "bandage/$filename"
                      else filename
                }

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_vg

    input:
    tuple val(sample), val(single_end), path(scaffolds), path (fasta) //from ch_metaspades_vg
    //path fasta from ch_fasta

    output:
    tuple val(sample), val(single_end), path("${sample}.vcf.gz*"), emit: ch_metaspades_vg_vcf
    path "*.bcftools_stats.txt", emit: ch_metaspades_vg_bcftools_mqc
    path "*.{gfa,png,svg}"

    script:
    """
    minimap2 -c -t $task.cpus -x asm20 $fasta $scaffolds > ${sample}.paf
    cat $scaffolds $fasta > ${sample}.withRef.fasta
    seqwish --paf-alns ${sample}.paf --seqs ${sample}.withRef.fasta --gfa ${sample}.gfa --threads $task.cpus
    vg view -Fv ${sample}.gfa --threads $task.cpus > ${sample}.vg
    vg convert -x ${sample}.vg > ${sample}.xg
    samtools faidx $fasta
    vg snarls ${sample}.xg > ${sample}.snarls
    for chrom in `cat ${fasta}.fai | cut -f1`
    do
        vg deconstruct -p \$chrom ${sample}.xg -r ${sample}.snarls --threads $task.cpus \\
            | bcftools sort -O v -T ./ \\
            | bgzip -c > ${sample}.\$chrom.vcf.gz
    done
    bcftools concat --output-type z --output ${sample}.vcf.gz *.vcf.gz
    tabix -p vcf -f ${sample}.vcf.gz
    bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
    if [ -s ${sample}.gfa ]
    then
        Bandage image ${sample}.gfa ${sample}.png --height 1000
        Bandage image ${sample}.gfa ${sample}.svg --height 1000
    fi
    """
}

process METASPADES_SNPEFF {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/metaspades/variants/snpeff", mode: params.publish_dir_mode

    //when:
    //!params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_vg && params.gff && !params.skip_snpeff

    input:
    tuple val(sample), val(single_end), path(vcf),file(db), file(config) //from ch_metaspades_vg_vcf
    //tuple file(db), file(config) //from ch_snpeff_db_metaspades
    val (index_base)

    output:
    path "*.snpEff.csv", emit: ch_metaspades_snpeff_mqc
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