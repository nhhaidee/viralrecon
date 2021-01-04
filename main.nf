#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/viralrecon
========================================================================================
 nf-core/viralrecon Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/viralrecon
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

include { helpMessage; checkHostname; nfcoreHeader;
          isOffline;
          validate_sample_sheet;
          CheckFasta} from './lib/helpers';

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

include { GUNZIP_FASTA; 
          GUNZIP_GFF; 
          UNTAR_KRAKEN2_DB; 
          CHECK_SAMPLESHEET;
          MAKE_BLAST_DB } from './modules/processes/preprocessing'
include { SRA_FASTQ_FTP; SRA_FASTQ_DUMP } from './modules/processes/downloadfiles'
include { CAT_FASTQ; FASTQC } from './modules/processes/fastqc'
include { FASTP } from './modules/processes/adaptertrimming'
include { BOWTIE2_INDEX; BOWTIE2 } from './modules/processes/bowtie'
include { SORT_BAM; SAMTOOLS_MPILEUP} from './modules/processes/samtools'
include { MAKE_SNPEFF_DB } from './modules/processes/snpeffdb'
include { get_software_versions; output_documentation } from './modules/processes/docs'
include { IVAR_TRIM; IVAR_VARIANTS; IVAR_CONSENSUS; IVAR_SNPEFF; IVAR_QUAST} from './modules/processes/ivar'
include { BCFTOOLS_VARIANTS; BCFTOOLS_CONSENSUS; BCFTOOLS_SNPEFF; BCFTOOLS_QUAST} from './modules/processes/bcf'
include { PICARD_MARKDUPLICATES; PICARD_METRICS } from './modules/processes/picard'
include { MOSDEPTH_GENOME; MOSDEPTH_AMPLICON; MOSDEPTH_AMPLICON_PLOT } from './modules/processes/mosdepth'
include { VARSCAN2; VARSCAN2_CONSENSUS; VARSCAN2_SNPEFF; VARSCAN2_QUAST } from './modules/processes/varscan'
include { BCFTOOLS_ISEC } from './modules/processes/intersect_ivar_bcf'
include { CUTADAPT } from './modules/processes/cutadapt'
include { KRAKEN2 } from './modules/processes/kraken2'
include { SPADES; SPADES_BLAST; SPADES_ABACAS } from './modules/processes/spades'
include { MULTIQC } from './modules/processes/multiqc'
include { FGBIO_TRIM_PRIMERS } from './modules/processes/primertrimming'


    // Has the run name been specified by the user?
    // this has the bonus effect of catching both -name and --name
workflow {

    custom_runName = params.name
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        custom_runName = workflow.runName
    }

    ////////////////////////////////////////////////////
    /* --          VALIDATE INPUTS                 -- */
    ////////////////////////////////////////////////////

    if (params.input) { sample_sheet_file = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

    if (params.trim_primer) { ch_primer_info = file(params.primer_info, checkIfExists: true) } else { exit 1, "Primer Infor file not specified!" }

    if (params.protocol != 'metagenomic' && params.protocol != 'amplicon') {
        exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'metagenomic' or 'amplicon'!"
    }

    if (params.protocol == 'amplicon' && !params.skip_assembly && !params.amplicon_fasta) {
        exit 1, "To perform de novo assembly in 'amplicon' mode please provide a valid amplicon fasta file!"
    }
    if (params.amplicon_fasta) { ch_amplicon_fasta = file(params.amplicon_fasta, checkIfExists: true) }

    if (params.protocol == 'amplicon' && !params.skip_variants && !params.amplicon_bed) {
        exit 1, "To perform variant calling in 'amplicon' mode please provide a valid amplicon BED file!"
    }
    if (params.amplicon_bed) { ch_amplicon_bed = file(params.amplicon_bed, checkIfExists: true) }

    callerList = [ 'varscan2', 'ivar', 'bcftools']
    callers = params.callers ? params.callers.split(',').collect{ it.trim().toLowerCase() } : []
    if ((callerList + callers).unique().size() != callerList.size()) {
        exit 1, "Invalid variant calller option: ${params.callers}. Valid options: ${callerList.join(', ')}"
    }

    assemblerList = [ 'spades', 'metaspades', 'unicycler', 'minia' ]
    assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
    if ((assemblerList + assemblers).unique().size() != assemblerList.size()) {
        exit 1, "Invalid assembler option: ${params.assemblers}. Valid options: ${assemblerList.join(', ')}"
    }

    // Viral reference files
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the Genome file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    }
    params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
    params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false

    if (params.fasta) {
        file(params.fasta, checkIfExists: true)
        
        lastPath = params.fasta.lastIndexOf(File.separator)
        lastExt = params.fasta.lastIndexOf(".")
        fasta_base = params.fasta.substring(lastPath+1)
        index_base = params.fasta.substring(lastPath+1,lastExt)
        if (params.fasta.endsWith('.gz')) {
            fasta_base = params.fasta.substring(lastPath+1,lastExt)
            index_base = fasta_base.substring(0,fasta_base.lastIndexOf("."))
        }
        
    } else {
        exit 1, "Viral genome fasta file not specified!"
    }

    ////////////////////////////////////////////////////
    /* --          CONFIG FILES                    -- */
    ////////////////////////////////////////////////////

    ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
    ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

    ////////////////////////////////////////////////////
    /* --          HEADER FILES                    -- */
    ////////////////////////////////////////////////////

    ch_blast_outfmt6_header = file("$baseDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
    ch_ivar_variants_header_mqc = file("$baseDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

    ////////////////////////////////////////////////////
    /* --                   AWS                    -- */
    ////////////////////////////////////////////////////

    // Check AWS batch settings
    if (workflow.profile.contains('awsbatch')) {
        // AWSBatch sanity checking
        if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
        // Check outdir paths to be S3 buckets if running on AWSBatch
        // related: https://github.com/nextflow-io/nextflow/issues/813
        if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
        // Prevent trace files to be stored on S3 since S3 does not support rolling files.
        if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                       HEADER LOG INFO                               -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    // Header log info
    log.info nfcoreHeader()
    def summary = [:]
    if (workflow.revision)               summary['Pipeline Release'] = workflow.revision
    summary['Run Name']                  = custom_runName ?: workflow.runName
    summary['Samplesheet']               = params.input
    summary['Protocol']                  = params.protocol
    if (params.protocol == 'amplicon')   summary['Amplicon Fasta File'] = params.amplicon_fasta
    if (params.protocol == 'amplicon')   summary['Amplicon BED File'] = params.amplicon_bed
    if (params.protocol == 'amplicon')   summary['Amplicon Left Suffix'] = params.amplicon_left_suffix
    if (params.protocol == 'amplicon')   summary['Amplicon Right Suffix'] = params.amplicon_right_suffix
    summary['Viral Genome']              = params.genome ?: 'Not supplied'
    summary['Viral Fasta File']          = params.fasta
    if (params.trim_primer)              summary['FLEX primer info '] = params.primer_info
    if (params.gff)                      summary['Viral GFF'] = params.gff
    if (params.save_reference)           summary['Save Genome Indices'] = 'Yes'
    if (params.save_sra_fastq)           summary['Save SRA FastQ'] = params.save_sra_fastq
    if (params.skip_sra)                 summary['Skip SRA Download'] = params.skip_sra
    if (!params.skip_adapter_trimming)  {
        if (params.cut_mean_quality)          summary['Fastp Mean Qual'] = params.cut_mean_quality
        if (params.qualified_quality_phred)   summary['Fastp Qual Phred'] = params.qualified_quality_phred
        if (params.unqualified_percent_limit) summary['Fastp Unqual % Limit'] = params.unqualified_percent_limit
        if (params.min_trim_length)           summary['Fastp Min Trim Length'] = params.min_trim_length
    } else {
        summary['Skip Adapter Trimming'] = 'Yes'
    }
    if (params.skip_amplicon_trimming)   summary['Skip Amplicon Trimming'] = 'Yes'
    if (params.save_trimmed)             summary['Save Trimmed'] = 'Yes'
    if (!params.skip_variants) {
        summary['Variant Calling Tools'] = params.callers
        summary['Min Mapped Reads']      = params.min_mapped_reads
        if (params.ivar_trim_noprimer)   summary['iVar Trim Exclude']  = 'Yes'
        summary['iVar Trim Min Len']     = params.ivar_trim_min_len
        summary['iVar Trim Min Qual']    = params.ivar_trim_min_qual
        summary['iVar Trim Window']      = params.ivar_trim_window_width
        if (params.filter_dups)          summary['Remove Duplicate Reads']  = 'Yes'
        if (params.filter_unmapped)      summary['Remove Unmapped Reads']  = 'Yes'
        summary['Mpileup Depth']         = params.mpileup_depth
        summary['Min Base Quality']      = params.min_base_qual
        summary['Min Read Depth']        = params.min_coverage
        summary['Min Allele Freq']       = params.min_allele_freq
        summary['Max Allele Freq']       = params.max_allele_freq
        if (params.varscan2_strand_filter) summary['Varscan2 Strand Filter'] = 'Yes'
        if (params.save_align_intermeds) summary['Save Align Intermeds'] =  'Yes'
        if (params.save_mpileup)         summary['Save mpileup'] = 'Yes'
        if (params.skip_markduplicates)  summary['Skip MarkDuplicates'] = 'Yes'
        if (params.skip_picard_metrics)  summary['Skip Picard Metrics'] = 'Yes'
        if (params.skip_mosdepth)        summary['Skip mosdepth'] = 'Yes'
        if (params.skip_snpeff)          summary['Skip SnpEff'] = 'Yes'
        if (params.skip_variants_quast)  summary['Skip Variants QUAST'] = 'Yes'
    } else {
        summary['Skip Variant Calling']  = 'Yes'
    }
    if (!params.skip_kraken2 && !params.skip_assembly) {
        if (params.kraken2_db)           summary['Host Kraken2 DB'] = params.kraken2_db
        if (params.kraken2_db_name)      summary['Host Kraken2 Name'] = params.kraken2_db_name
        if (params.kraken2_use_ftp)      summary['Kraken2 Use FTP'] = params.kraken2_use_ftp
        if (params.save_kraken2_fastq)   summary['Save Kraken2 FastQ'] = params.save_kraken2_fastq
    } else {
        summary['Skip Kraken2']          = 'Yes'
    }
    if (!params.skip_assembly) {
        summary['Assembly Tools']        = params.assemblers
        summary['Minia Kmer Size']       = params.minia_kmer
        if (params.skip_vg)              summary['Skip Variant Graph'] =  'Yes'
        if (params.skip_blast)           summary['Skip BLAST'] =  'Yes'
        if (params.skip_abacas)          summary['Skip ABACAS'] =  'Yes'
        if (params.skip_plasmidid)       summary['Skip PlasmidID'] =  'Yes'
        if (params.skip_assembly_quast)  summary['Skip Assembly QUAST'] =  'Yes'
    } else {
        summary['Skip Assembly']         = 'Yes'
    }
    if (params.skip_fastqc)              summary['Skip FastQC'] = 'Yes'
    if (params.skip_multiqc)             summary['Skip MultiQC'] = 'Yes'
    summary['Max Resources']             = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
    if (workflow.containerEngine)        summary['Container'] = "$workflow.containerEngine - $workflow.container"
    summary['Output dir']                = params.outdir
    summary['Publish dir mode']          = params.publish_dir_mode
    summary['Launch dir']                = workflow.launchDir
    summary['Working dir']               = workflow.workDir
    summary['Script dir']                = workflow.projectDir
    summary['User']                      = workflow.userName
    if (workflow.profile.contains('awsbatch')) {
        summary['AWS Region']            = params.awsregion
        summary['AWS Queue']             = params.awsqueue
        summary['AWS CLI']               = params.awscli
    }
    summary['Config Profile']            = workflow.profile
    if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
    if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
    if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
    if (params.email || params.email_on_fail) {
        summary['E-mail Address']        = params.email
        summary['E-mail on failure']     = params.email_on_fail
        summary['MultiQC maxsize']       = params.max_multiqc_email_size
    }
    log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
    log.info "-\033[2m--------------------------------------------------\033[0m-"

    // Check the hostnames against configured profiles
    checkHostname()


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                 UNZIP/UNTAR REFERENCE FILES                         -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    /*
    * PREPROCESSING: Uncompress genome fasta file
    */
    if (params.fasta) {
        fasta_file = file(params.fasta, checkIfExists: true)
        if (params.fasta.endsWith('.gz')) {
            ch_fasta = Channel.from(fasta_file) | GUNZIP_FASTA
        }
        else {
            ch_fasta = Channel.from(fasta_file)
        }
    }

    /*
    * PREPROCESSING: Uncompress gff annotation file
    */

    if (params.gff) {
        gff_file  = file(params.gff, checkIfExists: true)
        if (params.gff.endsWith('.gz')) {
            ch_gff = Channel.from(gff_file) | GUNZIP_GFF
        } 
        else {
            ch_gff = Channel.from(gff_file)
        }
    } 
    else {
        //See: https://nextflow-io.github.io/patterns/index.html#_optional_input
        ch_gff = file('NO_FILE')
    }

    /*
    * PREPROCESSING: Uncompress Kraken2 database
    */

    if (!params.skip_kraken2 && params.kraken2_db && !params.skip_assembly) {
        kraken2_db_file = file(params.kraken2_db, checkIfExists: true)
        if (params.kraken2_db.endsWith('.tar.gz')) {
            ch_kraken2_db = Channel.from(kraken2_db_file) | UNTAR_KRAKEN2_DB
        } 
        else {
            ch_kraken2_db = Channel.from(kraken2_db_file)
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                     PARSE DESIGN FILE                               -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


    /*
    * PREPROCESSING: Reformat samplesheet and check validity
    */

    ch_input = Channel.from(sample_sheet_file) \
            | CHECK_SAMPLESHEET \
            | splitCsv(header:true, sep:',')
            | map { validate_sample_sheet(it) }

    ch_reads_all = ch_input | map { it } 
    ch_reads_sra = ch_input | map { it }

    /*
    * STEP 1: Download and check SRA data
    */
    if (!params.skip_sra || !isOffline()) {

        ch_reads_sra_ftp = ch_reads_sra | filter { it[2] }
        ch_reads_sra_dump = ch_reads_sra | filter { it[2] }

        ch_sra_fastq_ftp = ch_reads_sra_ftp | SRA_FASTQ_FTP

        ch_sra_fastq_dump = ch_reads_sra_dump.map { it[0..3] } | SRA_FASTQ_DUMP

        ch_reads_all
            .filter { !it[2] }
            .concat(ch_sra_fastq_ftp, ch_sra_fastq_dump)
            .set { ch_reads_all }
    }

    ch_reads_all
        .map { [ it[0].split('_')[0..-2].join('_'), it[1], it[4] ] }
        .groupTuple(by: [0, 1])
        .map { [ it[0], it[1], it[2].flatten() ] }
        .set { ch_reads_all }

    /*
    * STEP 2: Merge FastQ files with the same sample identifier
    */

    ch_reads_all | CAT_FASTQ // emit ch_cat_fastqc and ch_cat_fastp

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                        FASTQ QC                                     -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    /*
    * STEP 3: FastQC on input reads after merging libraries from the same sample
    */

    ch_fastqc_raw_reports_mqc = CAT_FASTQ.out.ch_cat_fastq | FASTQC

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                        ADAPTER TRIMMING                             -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    /*
    * STEP 4: Fastp adapter trimming and quality filtering
    */
    if (!params.skip_adapter_trimming) {
        CAT_FASTQ.out.ch_cat_fastq | FASTP
    } 
    /*
    else {
        
        ch_cat_fastp
            .into { ch_fastp_bowtie2
                    ch_fastp_cutadapt
                    ch_fastp_kraken2 }
        ch_fastp_mqc = Channel.empty()
        ch_fastp_fastqc_mqc = Channel.empty()
    }
    */

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                  VARIANT CALLING PROCESSES                          -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    /*
    * PREPROCESSING: Build Bowtie2 index for viral genome
    */
    
    BOWTIE2_INDEX (ch_fasta, index_base)

    /*
    * PREPROCESSING: Build SnpEff database for viral genome
    */
    if ((!params.skip_variants || !params.skip_assembly) && params.gff  && !params.skip_snpeff) {
        MAKE_SNPEFF_DB(ch_fasta, ch_gff, index_base)
    }
        

    /*
    * STEP 5.1: Map read(s) with Bowtie 2
    */
    FASTP.out.ch_fastp.combine(BOWTIE2_INDEX.out.ch_index)
    if (!params.skip_adapter_trimming){ 
        ch_bowtie2_align = FASTP.out.ch_fastp.combine(BOWTIE2_INDEX.out.ch_index)
        BOWTIE2 (ch_bowtie2_align, index_base)
    }
    else{
        ch_bowtie2_align = CAT_FASTQ.out.ch_cat_fastq.combine(BOWTIE2_INDEX.out.ch_index)
        BOWTIE2 (CAT_FASTQ.out.ch_cat_fastq, index_base)
    }

    /*
    Add Primer Trimming Step using FGBIO
    */
    if (params.trim_primer){

        /*
        Perform FGBIO Trim Primer
        */
        FGBIO_TRIM_PRIMERS(ch_primer_info, BOWTIE2.out.ch_bowtie2_bam)
    }
    
    /*
    * STEP 5.2: Convert BAM to coordinate sorted BAM
    */
    if (params.trim_primer){
        SORT_BAM(FGBIO_TRIM_PRIMERS.out.ch_fgbio_trim_primers_bam)
    }
    else{
        SORT_BAM(BOWTIE2.out.ch_bowtie2_bam)
    }


    //////// Remove samples that failed mapped read threshold (later) ///////


    /*
    * STEP 5.3: Trim amplicon sequences with iVar
    */
    if (params.protocol == 'amplicon'){
        ch_sort_bam = SORT_BAM.out.ch_sort_bam | map {it[0..2]}
        IVAR_TRIM(ch_sort_bam, ch_amplicon_bed)
    }

    /*
    * STEP 5.4: Picard MarkDuplicates
    */
    ch_picard_markduplicates = IVAR_TRIM.out.ch_ivar_trim_bam.combine(ch_fasta)
    ch_picard_markduplicates | PICARD_MARKDUPLICATES

    /*
    * STEP 5.5: Picard CollectMultipleMetrics and CollectWgsMetrics
    */
    ch_picard_metrics = PICARD_MARKDUPLICATES.out.ch_markdup_bam.combine(ch_fasta)
    ch_picard_metrics | PICARD_METRICS

    /*
    * STEP 5.6.1: mosdepth genome-wide coverage
    */
    PICARD_MARKDUPLICATES.out.ch_markdup_bam | MOSDEPTH_GENOME

    /*
    * STEP 5.6.2: mosdepth amplicon coverage and plots
    */

    if (params.protocol == 'amplicon') {

        MOSDEPTH_AMPLICON(PICARD_MARKDUPLICATES.out.ch_markdup_bam, ch_amplicon_bed)

        MOSDEPTH_AMPLICON.out.ch_mosdepth_amplicon_region_bed.collect() | MOSDEPTH_AMPLICON_PLOT
    }

    ////////////////////////////////////////////////////
    /* --              VARSCAN2                    -- */
    ////////////////////////////////////////////////////

    /*
    * STEP 5.7: Create mpileup file for all variant callers
    */
    ch_samtools_mpileup  = PICARD_MARKDUPLICATES.out.ch_markdup_bam.combine(ch_fasta)
    ch_samtools_mpileup | SAMTOOLS_MPILEUP

    /*
    * STEP 5.7.1: Variant calling with VarScan 2
    */
    ch_varscan2  = SAMTOOLS_MPILEUP.out.ch_mpileup.combine(ch_fasta)
    ch_varscan2 | VARSCAN2


    if (!params.skip_variants && 'varscan2' in callers){

        /*
        * STEP 5.7.1.1: Genome consensus generation with BCFtools and masked with BEDTools
        */
        ch_varscan2_consensus_input = PICARD_MARKDUPLICATES.out.ch_markdup_bam.join(VARSCAN2.out.ch_varscan2_highfreq, by: [0,1]).combine(ch_fasta)
        ch_varscan2_consensus_input | VARSCAN2_CONSENSUS

        /*
        * STEP 5.7.1.2: VarScan 2 variant calling annotation with SnpEff and SnpSift
        */

        if(params.gff && !params.skip_snpeff) {

            ch_varscan2_snpeff_input = VARSCAN2.out.ch_varscan2_highfreq.join(VARSCAN2.out.ch_varscan2_lowfreq, by: [0,1]).combine(MAKE_SNPEFF_DB.out.ch_snpeff_db)
            VARSCAN2_SNPEFF(ch_varscan2_snpeff_input, index_base)

        }

        /*
        * STEP 5.7.1.3: VarScan 2 consensus sequence report with QUAST
        */
        if (!params.skip_variants_quast){

            VARSCAN2_QUAST(VARSCAN2_CONSENSUS.out.ch_varscan2_consensus.collect{ it[2]}, ch_fasta, ch_gff)

        }
    }

    ////////////////////////////////////////////////////
    /* --                IVAR                      -- */
    ////////////////////////////////////////////////////

    /*
    * STEP 5.7.2: Variant calling with iVar
    */
    if (!params.skip_variants && 'ivar' in callers){

        ch_ivar_variants_input = SAMTOOLS_MPILEUP.out.ch_mpileup.combine(ch_fasta).combine(ch_gff)
        IVAR_VARIANTS(ch_ivar_variants_input, ch_ivar_variants_header_mqc)
        
        /*
        * STEP 5.7.2.1: Generate consensus sequence with iVar
        */
        ch_ivar_consensus_input = SAMTOOLS_MPILEUP.out.ch_mpileup.combine(ch_fasta)
        IVAR_CONSENSUS(ch_ivar_consensus_input)


        /*
        * STEP 5.7.2.2: iVar variant calling annotation with SnpEff and SnpSift
        */
        if (params.gff && !params.skip_snpeff)
        {
            ch_ivar_snpeff_input = IVAR_VARIANTS.out.ch_ivar_highfreq.join(IVAR_VARIANTS.out.ch_ivar_lowfreq, by: [0,1]).combine(MAKE_SNPEFF_DB.out.ch_snpeff_db)
            IVAR_SNPEFF(ch_ivar_snpeff_input, index_base)
        }
        /*
        * STEP 5.7.2.3: iVar consensus sequence report with QUAST
        */
        if(!params.skip_variants_quast){

            IVAR_QUAST(IVAR_CONSENSUS.out.ch_ivar_consensus.collect{ it[2]}, ch_fasta, ch_gff)
        }

    }

    ////////////////////////////////////////////////////
    /* --              BCFTOOLS                    -- */
    ////////////////////////////////////////////////////


    if (!params.skip_variants && 'bcftools' in callers){

        /*
        * STEP 5.7.3: Variant calling with BCFTools
        */
        BCFTOOLS_VARIANTS(PICARD_MARKDUPLICATES.out.ch_markdup_bam.combine(ch_fasta))

        /*
        * STEP 5.7.3.1: Genome consensus generation with BCFtools and masked with BEDTools
        */
        BCFTOOLS_CONSENSUS(PICARD_MARKDUPLICATES.out.ch_markdup_bam.join(BCFTOOLS_VARIANTS.out.ch_bcftools_variants, by: [0,1]).combine(ch_fasta), index_base)

        /*
        * STEP 5.7.3.2: BCFTools variant calling annotation with SnpEff and SnpSift
        */
        if (params.gff && !params.skip_snpeff){
           BCFTOOLS_SNPEFF(BCFTOOLS_VARIANTS.out.ch_bcftools_variants.combine(MAKE_SNPEFF_DB.out.ch_snpeff_db), index_base) 
        }
        /*
        * STEP 5.7.3.3: BCFTools consensus sequence report with QUAST
        */
        if (!params.skip_variants_quast){
            BCFTOOLS_QUAST(BCFTOOLS_CONSENSUS.out.ch_bcftools_consensus_masked.collect{ it[2] }, ch_fasta, ch_gff)
        }

    }
    ////////////////////////////////////////////////////
    /* --            INTERSECT VARIANTS            -- */
    ////////////////////////////////////////////////////

    /*
    * STEP 5.8: Intersect variants with BCFTools
    */
    if (!params.skip_variants && callers.size() > 2) {

        ch_intersect_ivar_bcf = VARSCAN2.out.ch_varscan2_highfreq.join(IVAR_VARIANTS.out.ch_ivar_highfreq, by: [0,1]).join(BCFTOOLS_VARIANTS.out.ch_bcftools_variants, by: [0,1])

        ch_intersect_ivar_bcf | BCFTOOLS_ISEC
    }
   

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /* --                                                                     -- */
    /* --                          MULTIQC                                    -- */
    /* --                                                                     -- */
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    Channel.from(summary.collect{ [it.key, it.value] })
        .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
        .reduce { a, b -> return [a, b].join("\n            ") }
        .map { x -> """
        id: 'nf-core-viralrecon-summary'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/viralrecon Workflow Summary'
        section_href: 'https://github.com/nf-core/viralrecon'
        plot_type: 'html'
        data: |
            <dl class=\"dl-horizontal\">
                $x
            </dl>
        """.stripIndent() }
        .set { ch_workflow_summary }

    
    get_software_versions()
    if (!params.skip_multiqc) {
        MULTIQC(
            custom_runName,
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            FASTQC.out.ch_fastqc_raw_reports_mqc.collect().ifEmpty([]),
            FASTP.out.ch_fastp_mqc.collect().ifEmpty([]),
            FASTP.out.ch_fastp_fastqc_mqc.collect().ifEmpty([]),
            BOWTIE2.out.ch_bowtie2_mqc.collect().ifEmpty([]),
            SORT_BAM.out.ch_sort_bam_flagstat_mqc.collect().ifEmpty([]),
            IVAR_TRIM.out.ch_ivar_trim_flagstat_mqc.collect().ifEmpty([]),
            IVAR_TRIM.out.ch_ivar_trim_log_mqc.collect().ifEmpty([]),
            PICARD_MARKDUPLICATES.out.ch_markdup_bam_flagstat_mqc.collect().ifEmpty([]),
            PICARD_MARKDUPLICATES.out.ch_markdup_bam_metrics_mqc.collect().ifEmpty([]),
            PICARD_METRICS.out.ch_picard_metrics_mqc.collect().ifEmpty([]),
            MOSDEPTH_GENOME.out.ch_mosdepth_genome_mqc.collect().ifEmpty([]),
            VARSCAN2.out.ch_varscan2_log_mqc.collect().ifEmpty([]),
            VARSCAN2.out.ch_varscan2_bcftools_highfreq_mqc.collect().ifEmpty([]),
            VARSCAN2_SNPEFF.out.ch_varscan2_snpeff_highfreq_mqc.collect().ifEmpty([]),
            VARSCAN2_QUAST.out.ch_varscan2_quast_mqc.collect().ifEmpty([]),
            IVAR_VARIANTS.out.ch_ivar_count_mqc.collect().ifEmpty([]),
            IVAR_VARIANTS.out.ch_ivar_bcftools_highfreq_mqc.collect().ifEmpty([]),
            IVAR_SNPEFF.out.ch_ivar_snpeff_highfreq_mqc.collect().ifEmpty([]),
            IVAR_QUAST.out.ch_ivar_quast_mqc.collect().ifEmpty([]),
            BCFTOOLS_VARIANTS.out.ch_bcftools_variants_mqc.collect().ifEmpty([]),
            BCFTOOLS_SNPEFF.out.ch_bcftools_snpeff_mqc.collect().ifEmpty([]),
            BCFTOOLS_QUAST.out.ch_bcftools_quast_mqc.collect().ifEmpty([])
        )
    }
    output_documentation(ch_output_docs,ch_output_docs_images)
}
