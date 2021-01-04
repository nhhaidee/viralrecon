
process BOWTIE2_INDEX {
    tag "$fasta"
    label 'process_medium'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    when:
    !params.skip_variants

    input:
    path fasta /* from ch_fasta */
    val (index_base)

    output:
    path "Bowtie2Index", emit: ch_index

    script:    
    """
    bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $fasta \\
        $index_base
    mkdir Bowtie2Index && mv ${index_base}* Bowtie2Index
    """
}

process BOWTIE2 {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".log")) "log/$filename"
                      else params.save_align_intermeds ? filename : null
                }

    when:
    !params.skip_variants

    input:
    tuple val(sample), val(single_end), path(reads), path(index) /* from ch_fastp_bowtie2*/
     /* from ch_index*/
    val (index_base)

    output:
    tuple val(sample), val(single_end), path("*.bam"), emit: ch_bowtie2_bam /* into ch_bowtie2_bam */
    path "*.log", emit: ch_bowtie2_mqc /*into ch_bowtie2_mqc*/

    script:
    input_reads = single_end ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    filter = params.filter_unmapped ? "-F4" : ""
    """
    bowtie2 \\
        --threads $task.cpus \\
        --local \\
        --very-sensitive-local \\
        -x ${index}/${index_base} \\
        $input_reads \\
        2> ${sample}.bowtie2.log \\
        | samtools view -@ $task.cpus -b -h -O BAM -o ${sample}.bam $filter -
    """
}