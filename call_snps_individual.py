#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=8192
import sys, os
from ast import literal_eval
from utils import *
from utils.configuration import *
import tempfile
import glob


#====================#
# Load Configuration #
#====================#

bam = dotdictify(literal_eval(sys.argv[2]))
cf = config(sys.argv[1])

#==========#
# BCFTools #
#==========#
if 'bcftools' in cf.snps:
    print "GREAT"
    # Test to see if a union set of variants has been constructed
    call_individual = file_exists(cf.union_variants["bcftools"]["ALL"])
    if call_individual:
        vcf_filename = bam.vcf_files.bcftools_union
        vcf_form = "union"
        union_variant_file = " -T {cf.union_variants.bcftools.ALL} ".format(**locals())
        output_all_sites = True
    else:
        vcf_filename = bam.vcf_files.bcftools_individual
        vcf_form = "individual"
        union_variant_file = ""
        output_all_sites = False

    print call_individual, "CALL INDIVIDUAL"

    samtools_mpileup = "samtools mpileup {cf.snps.samtools_options} -g -f {cf.reference_file} -r __region__ {bam.bam_merged_filename}".format(**locals())
    bcftools_call = "bcftools call {cf.snps.bcftools_options} {union_variant_file} ".format(**locals())
    # When calling union vcf files, output all sites, not just variant ones.
    if output_all_sites == True:
        bcftools_call = bcftools_call.replace(" -v ", "")
    xarg_command = "REGION='__region__'; {samtools_mpileup} | {bcftools_call} -O z > {cf.vcf_dir}/{bam.SM}.TMP.${{REGION/:/_}}.bcftools.{vcf_form}.vcf.gz && bcftools index {cf.vcf_dir}/{bam.SM}.TMP.${{REGION/:/_}}.bcftools.{vcf_form}.vcf.gz".format(**locals())
    # Replace last colon (screws up file names)
    chrom_chunks = cf.chunk_genome()
    chrom_chunks_joined = ','.join(chrom_chunks)
    bcftools = """{echo} -n {chrom_chunks_joined} | {xargs}  -t -d ',' -P 24 -I {{}} bash -c '{xarg_command}' """.format(**locals()).replace("__region__","{}")
    cf.command(bcftools)

    # Split out snps and indels
    split_indels = "bcftools view --include='INDEL=0' {union_variant_file}"

    #=============#
    # SNP Filters #
    #=============#
    # Heterozygous Polarization
    if cf.snps.bcftools.__heterozygous_polarization == True:
        snp_filters = "| bcftools view -O v | python {script_dir}/het_polarization.py ".format(**locals())
    else:
        snp_filters = ""

    soft_filters = cf.snps.bcftools.__soft_filters
    hard_filters = cf.snps.bcftools.__hard_filters
    snp_filters += construct_filters(cf.snps.bcftools.__soft_filters)
    snp_filters += construct_filters(cf.snps.bcftools.__hard_filters, soft=False)

    #===============#
    # Indel Filters #
    #===============#

    if cf.indels.bcftools.__heterozygous_polarization == True:
        indel_filters = "| bcftools view -O v | python {script_dir}/het_polarization.py ".format(**locals())
    else:
        indel_filters = ""

    soft_filters = cf.snps.bcftools.__soft_filters
    hard_filters = cf.snps.bcftools.__hard_filters
    indel_filters += construct_filters(cf.snps.bcftools.__soft_filters)
    indel_filters += construct_filters(cf.indels.bcftools.__hard_filters, soft=False)

    # Merge, sort, and index
    concat_files = ["{cf.vcf_dir}/{bam.SM}.TMP.{chunk}.bcftools.{vcf_form}.vcf.gz".format(**locals()).replace(":","_") for chunk in chrom_chunks]
    concat_list = ' '.join(concat_files)

    # Put together SNP list

    bcftools_concat = """bcftools concat -O z {concat_list}  > {vcf_filename.ALL};
                         bcftools index -f {vcf_filename.ALL};
                         bcftools view -O z {vcf_filename.ALL} --include='INDEL=0' {snp_filters} > {vcf_filename.SNP};
                         bcftools index -f {vcf_filename.SNP};
                         bcftools view -O z {vcf_filename.ALL} --include='INDEL=1' {indel_filters} > {vcf_filename.INDEL};
                         bcftools index -f {vcf_filename.INDEL}""".format(**locals()) # Add indel filters.

    cf.command(bcftools_concat)

    # Remove temporary files
    if check_seq_file(vcf_filename.SNP) and file_exists(vcf_filename.INDEL):
        comm = "rm " + ' '.join(concat_files + [x + ".csi" for x in concat_files])
        cf.command(comm)

