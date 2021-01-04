process SPADES {
    tag "$sample"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/spades", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".png")) "bandage/$filename"
                      else if (filename.endsWith(".svg")) "bandage/$filename"
                      else filename
                }

    //when:
    //!params.skip_assembly && 'spades' in assemblers

    input:
    tuple val(sample), val(single_end), path(reads)/* from ch_kraken2_spades*/

    output:
    tuple val(sample), val(single_end), path("*scaffolds.fa"), emit: ch_spades /*into ch_spades_blast,
                                                                   ch_spades_abacas,
                                                                   ch_spades_plasmidid,
                                                                   ch_spades_quast,
                                                                   ch_spades_vg*/
    path "*assembly.{gfa,png,svg}"


    script:
    input_reads = single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    spades.py \\
        --threads $task.cpus \\
        $input_reads \\
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

process SPADES_BLAST {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/spades/blast", mode: params.publish_dir_mode

    //when:
    //!params.skip_assembly && 'spades' in assemblers && !params.skip_blast

    input:
    tuple val(sample), val(single_end), path(scaffold), path (db)/* from ch_spades_blast */
    /* path db rom ch_blast_db */
    path header /*from ch_blast_outfmt6_header*/
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

process SPADES_ABACAS {
    tag "$sample"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/spades/abacas", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
                      else filename
                }

    //when:
    //!params.skip_assembly && 'spades' in assemblers && !params.skip_abacas

    input:
    tuple val(sample), val(single_end), path(scaffold), path(fasta) /*from ch_spades_abacas */
    /*path fasta from ch_fasta*/

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