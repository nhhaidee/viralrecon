#!/usr/bin/env python

import os
import sys
import errno
import argparse
import yaml


def parse_args(args=None):
    Description = 'Create custom spreadsheet for pertinent MultiQC metrics generated by the nf-core/viralrecon pipeline.'
    Epilog = "Example usage: python multiqc_to_custom_tsv.py"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-md', '--multiqc_data_dir', type=str, dest="MULTIQC_DATA_DIR", default='multiqc_data', help="Full path to directory containing YAML files for each module, as generated by MultiQC. (default: 'multiqc_data').")
    parser.add_argument('-of', '--out_file', type=str, dest="OUT_FILE", default='summary_metrics.tsv', help="Full path to output file (default: 'summary_metrics.tsv').")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


# Find key in dictionary created from YAML file recursively
# From https://stackoverflow.com/a/37626981
def find_tag(d, tag):
    if tag in d:
        yield d[tag]
    for k,v in d.items():
        if isinstance(v, dict):
            for i in find_tag(v, tag):
                yield i


def yaml_fields_to_dict(YAMLFile,AppendDict={},FieldMappingList=[]):
    with open(YAMLFile) as f:
        yaml_dict = yaml.safe_load(f)
        for k in yaml_dict.keys():
            key = k
            if YAMLFile.find('multiqc_picard_insertSize') != -1:
                key = k[:-3]
            if key not in AppendDict:
                AppendDict[key] = {}
            if FieldMappingList != []:
                for i,j in FieldMappingList:
                    val = list(find_tag(yaml_dict[k], j[0]))
                    if len(val) != 0:
                        val = val[0]
                        if len(j) == 2:
                            val = list(find_tag(val, j[1]))[0]
                        if j[0] in ['number_of_SNPs', 'number_of_indels', 'MISSENSE']:
                            val = int(val)
                        if i not in AppendDict[key]:
                            AppendDict[key][i] = val
                        else:
                            print('WARNING: {} key already exists in dictionary so will be overwritten. YAML file {}.'.format(i,YAMLFile))
            else:
                AppendDict[key] = yaml_dict[k]
    return AppendDict


def metrics_dict_to_file(FileFieldList,MultiQCDataDir,OutFile,SectionHeader='## METRICS',append=False):
    MetricsDict = {}
    FieldList = []
    for yamlFile,mappingList in FileFieldList:
        yamlFile = os.path.join(MultiQCDataDir,yamlFile)
        if os.path.exists(yamlFile):
            MetricsDict = yaml_fields_to_dict(YAMLFile=yamlFile,AppendDict=MetricsDict,FieldMappingList=mappingList)
            FieldList += [x[0] for x in mappingList]
    if MetricsDict != {}:
        make_dir(os.path.dirname(OutFile))
        if append:
            fout = open(OutFile,'a')
        else:
            fout = open(OutFile,'w')
        header = ['Sample'] + FieldList
        fout.write('{}\n'.format(SectionHeader))
        fout.write('{}\n'.format('\t'.join(header)))
        for k in sorted(MetricsDict.keys()):
            rowList = [k]
            for field in FieldList:
                if field in MetricsDict[k]:
                    rowList.append(MetricsDict[k][field])
                else:
                    rowList.append('NA')
            fout.write('{}\n'.format('\t'.join(map(str,rowList))))
        fout.close()
    return MetricsDict


def main(args=None):
    args = parse_args(args)

    ## File names for MultiQC YAML along with fields to fetch from each file
    VariantFileFieldList = [
        ('multiqc_fastp.yaml',                                     [('# Input reads', ['before_filtering','total_reads']),
                                                                    ('# Trimmed reads (fastp)', ['after_filtering','total_reads'])]),
        ('multiqc_samtools_flagstat_samtools_bowtie2.yaml',        [('% Mapped reads (viral)', ['mapped_passed_pct'])]),
        ('multiqc_ivar_summary.yaml',                              [('# Reads trimmed (iVar)', ['trimmed_reads'])]),
        ('multiqc_samtools_flagstat_samtools_ivar.yaml',           [('# Trimmed reads (iVar)', ['flagstat_total'])]),
        ('multiqc_samtools_flagstat_samtools_markduplicates.yaml', [('# Duplicate reads', ['duplicates_passed']),
                                                                    ('# Reads after MarkDuplicates', ['flagstat_total'])]),
        ('multiqc_picard_insertSize.yaml',                         [('Insert size mean', ['MEAN_INSERT_SIZE']),
                                                                    ('Insert size std dev', ['STANDARD_DEVIATION'])]),
        ('multiqc_picard_wgsmetrics.yaml',                         [('Coverage mean', ['MEAN_COVERAGE']),
                                                                    ('Coverage std dev', ['SD_COVERAGE']),
                                                                    ('% Coverage > 10x', ['PCT_10X'])]),
        ('multiqc_bcftools_stats_bcftools_varscan2.yaml',          [('# High conf SNPs (VarScan 2)', ['number_of_SNPs']),
                                                                    ('# High conf INDELs (VarScan 2)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_varscan2.yaml',                    [('# Missense variants (VarScan 2)', ['MISSENSE'])]),
        ('multiqc_quast_quast_varscan2.yaml',                      [('# Ns per 100kb consensus (VarScan 2)', ["# N's per 100 kbp"])]),
        ('multiqc_bcftools_stats_bcftools_ivar.yaml',              [('# High conf SNPs (iVar)', ['number_of_SNPs']),
                                                                    ('# High conf INDELs (iVar)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_ivar.yaml',                        [('# Missense variants (iVar)', ['MISSENSE'])]),
        ('multiqc_quast_quast_ivar.yaml',                          [('# Ns per 100kb consensus (iVar)', ["# N's per 100 kbp"])]),
        ('multiqc_bcftools_stats_bcftools_bcftools.yaml',          [('# High conf SNPs (BCFTools)', ['number_of_SNPs']),
                                                                    ('# High conf INDELs (BCFTools)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_bcftools.yaml',                    [('# Missense variants (BCFTools)', ['MISSENSE'])]),
        ('multiqc_quast_quast_bcftools.yaml',                      [('# Ns per 100kb consensus (BCFTools)', ["# N's per 100 kbp"])]),
    ]

    AssemblyFileFieldList = [
        ('multiqc_fastp.yaml',                                     [('# Input reads', ['before_filtering','total_reads'])]),
        #('multiqc_fastqc_fastqc_cutadapt.yaml',                    [('# Trimmed reads (Cutadapt)', ['Total Sequences'])]), ## HAVE TO MULTIPLY BY 2 FOR PE READS?
        ('multiqc_quast_quast_spades.yaml',                        [('# Contigs (SPAdes)', ['# contigs']),
                                                                    ('# Contigs > 5kb (SPAdes)', ['# contigs (>= 5000 bp)']),
                                                                    ('Largest contig (SPAdes)', ['Largest contig']),
                                                                    ('% Genome fraction (SPAdes)', ['Genome fraction (%)']),
                                                                    ('N50 (SPAdes)', ['N50'])]),
        ('multiqc_bcftools_stats_bcftools_spades.yaml',            [('# SNPs (SPAdes)', ['number_of_SNPs']),
                                                                    ('# INDELs (SPAdes)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_spades.yaml',                      [('# Missense variants (SPAdes)', ['MISSENSE'])]),
        ('multiqc_quast_quast_metaspades.yaml',                    [('# Contigs (metaSPAdes)', ['# contigs']),
                                                                    ('# Contigs > 5kb (metaSPAdes)', ['# contigs (>= 5000 bp)']),
                                                                    ('Largest contig (metaSPAdes)', ['Largest contig']),
                                                                    ('% Genome fraction (metaSPAdes)', ['Genome fraction (%)']),
                                                                    ('N50 (metaSPAdes)', ['N50'])]),
        ('multiqc_bcftools_stats_bcftools_metaspades.yaml',        [('# SNPs (metaSPAdes)', ['number_of_SNPs']),
                                                                    ('# INDELs (metaSPAdes)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_metaspades.yaml',                  [('# Missense variants (metaSPAdes)', ['MISSENSE'])]),
        ('multiqc_quast_quast_unicycler.yaml',                     [('# Contigs (Unicycler)', ['# contigs']),
                                                                    ('# Contigs > 5kb (Unicycler)', ['# contigs (>= 5000 bp)']),
                                                                    ('Largest contig (Unicycler)', ['Largest contig']),
                                                                    ('% Genome fraction (Unicycler)', ['Genome fraction (%)']),
                                                                    ('N50 (Unicycler)', ['N50'])]),
        ('multiqc_bcftools_stats_bcftools_unicycler.yaml',         [('# SNPs (Unicycler)', ['number_of_SNPs']),
                                                                    ('# INDELs (Unicycler)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_unicycler.yaml',                   [('# Missense variants (Unicycler)', ['MISSENSE'])]),
        ('multiqc_quast_quast_minia.yaml',                         [('# Contigs (minia)', ['# contigs']),
                                                                    ('# Contigs > 5kb (minia)', ['# contigs (>= 5000 bp)']),
                                                                    ('Largest contig (minia)', ['Largest contig']),
                                                                    ('% Genome fraction (minia)', ['Genome fraction (%)']),
                                                                    ('N50 (minia)', ['N50'])]),
        ('multiqc_bcftools_stats_bcftools_minia.yaml',             [('# SNPs (minia)', ['number_of_SNPs']),
                                                                    ('# INDELs (minia)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_minia.yaml',                       [('# Missense variants (minia)', ['MISSENSE'])])
    ]

    ## Write variant calling metrics to file
    MetricsDict = metrics_dict_to_file(FileFieldList=VariantFileFieldList,
                                       MultiQCDataDir=args.MULTIQC_DATA_DIR,
                                       OutFile=args.OUT_FILE,
                                       SectionHeader='## METRICS: VARIANT CALLING',
                                       append=False)

    ## Write de novo assembly metrics to file
    append = True
    SectionHeader = '\n## METRICS: DE NOVO ASSEMBLY'
    if MetricsDict == {}:
        append = False
        SectionHeader = SectionHeader.strip()
    metrics_dict_to_file(FileFieldList=AssemblyFileFieldList,
                         MultiQCDataDir=args.MULTIQC_DATA_DIR,
                         OutFile=args.OUT_FILE,
                         SectionHeader=SectionHeader,
                         append=append)


if __name__ == '__main__':
    sys.exit(main())
