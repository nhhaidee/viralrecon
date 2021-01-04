def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

process GUNZIP_FASTA {
    label 'error_retry'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    input:
    path fasta

    output:
    path "$unzip"

    script:
    unzip = fasta.toString() - '.gz'
    """
    pigz -f -d -p $task.cpus $fasta
    """
}

process GUNZIP_GFF {
    label 'error_retry'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    input:
    path gff

    output:
    path "$unzip"

    script:
    unzip = gff.toString() - '.gz'
    """
    pigz -f -d -p $task.cpus $gff
    """
}

process UNTAR_KRAKEN2_DB {
    label 'error_retry'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    input:
    path db

    output:
    path "$untar"

    script:
    untar = db.toString() - '.tar.gz'
    """
    tar -xvf $db
    """
}

process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "preprocess/sra/$filename"
                      else "pipeline_info/$filename"
                }

    input:
    path samplesheet

    output:
    path "samplesheet.valid.csv"/* into ch_samplesheet_reformat, emit: samplesheet_results*/
    //path "sra_run_info.tsv" optional true

    script:  // These scripts are bundled with the pipeline, in nf-core/viralrecon/bin/
    run_sra = !params.skip_sra && !isOffline()
    """
    awk -F, '{if(\$1 != "" && \$2 != "") {print \$0}}' $samplesheet > nonsra_id.csv
    check_samplesheet.py nonsra_id.csv nonsra.samplesheet.csv

    awk -F, '{if(\$1 != "" && \$2 == "" && \$3 == "") {print \$1}}' $samplesheet > sra_id.list
    if $run_sra && [ -s sra_id.list ]
    then
        fetch_sra_runinfo.py sra_id.list sra_run_info.tsv --platform ILLUMINA --library_layout SINGLE,PAIRED
        sra_runinfo_to_samplesheet.py sra_run_info.tsv sra.samplesheet.csv
    fi

    if [ -f nonsra.samplesheet.csv ]
    then
        head -n 1 nonsra.samplesheet.csv > samplesheet.valid.csv
    else
        head -n 1 sra.samplesheet.csv > samplesheet.valid.csv
    fi
    tail -n +2 -q *sra.samplesheet.csv >> samplesheet.valid.csv
    """
}

process MAKE_BLAST_DB {
    tag "$fasta"
    label 'process_medium'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    //when:
    //!params.skip_assembly && !params.skip_blast

    input:
    path fasta /*from ch_fasta */

    output:
    path "BlastDB", emit: ch_blast_db

    script:
    """
    makeblastdb \\
        -in $fasta \\
        -parse_seqids \\
        -dbtype nucl
    mkdir BlastDB && mv ${fasta}* BlastDB
    """
}

process KRAKEN2_BUILD {
    tag "$db"
    label 'process_high'
    if (params.save_reference) {
        publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
    }

    //when:
    //!params.skip_assembly

    output:
    path "$db", emit: ch_kraken2_db

    script:
    db = "kraken2_${params.kraken2_db_name}"
    ftp = params.kraken2_use_ftp ? "--use-ftp" : ""
    """
    kraken2-build --db $db --threads $task.cpus $ftp --download-taxonomy
    kraken2-build --db $db --threads $task.cpus $ftp --download-library $params.kraken2_db_name
    kraken2-build --db $db --threads $task.cpus $ftp --build
    """
}

