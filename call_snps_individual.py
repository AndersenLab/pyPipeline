#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
import tempfile
import glob


#====================#
# Load Configuration #
#====================#

bam = sys.argv[2]
SM = bam.replace(".bam","")
config, log, c_log = load_config_and_log(sys.argv[1], "snps")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
snps = COMMANDS.snps # Pulls out snp types.
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]

vcf_dir = "{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals())
makedir("{vcf_dir}".format(**locals()))


# Setup Chromosome Chunks
chrom_chunks = [x for x in chunk_genome(OPTIONS.chrom_chunk_kb*1000, reference)]
chrom_list = '\n'.join([x for x in chunk_genome(OPTIONS.chrom_chunk_kb*1000, reference)])
chrom_chunks_file = "{vcf_dir}_chrom_chunks.txt".format(**locals())
with open(chrom_chunks_file, "w+") as f:
	f.write(chrom_list)
	
#==========#
# BCFTools #
#==========#
if 'bcftools' in snps:
	# Test to see if a union set of variants has been constructed
	ind_union_variant_set = "{OPTIONS.analysis_dir}/bcftools.{OPTIONS.union_variants}.txt".format(**locals())

	complete_individual = "{vcf_dir}/{SM}.bcftools.individual.vcf.gz".format(**locals())
	individual_vcfs_check = (not file_exists(complete_individual) and COMMANDS.snps.snp_options.remove_temp == True)
	if individual_vcfs_check and file_exists(ind_union_variant_set):
		ind_union_filename = "union"
		call_union_sites = " -T {ind_union_variant_set} ".format(**locals())
		output_all_sites = True
	else:
		ind_union_filename = "individual"
		call_union_sites = ""
		output_all_sites = False

	# Completed File Name
	complete_ind_vcf = "{vcf_dir}/{SM}.bcftools.{ind_union_filename}.vcf.gz".format(**locals())
	if not file_exists(complete_ind_vcf):
		samtools_options = format_command(snps["samtools"])
		bcftools_options = format_command(snps["bcftools"])
		samtools_mpileup = "samtools mpileup {samtools_options} -g -f {reference} -r __region__ {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{bam}".format(**locals())
		
		bcftools_call = "bcftools call --skip-variants indels {bcftools_options} {call_union_sites} ".format(**locals())

		# When calling union vcf files, output all sites, not just variant ones.
		if output_all_sites == True:
			bcftools_call = bcftools_call.replace(" -v ", "")

		xarg_command = "REGION='__region__'; {samtools_mpileup} | {bcftools_call} -O z > {vcf_dir}/TMP.{SM}.${{REGION/:/_}}.bcftools.{ind_union_filename}.vcf.gz && bcftools index {vcf_dir}/TMP.{SM}.${{REGION/:/_}}.bcftools.{ind_union_filename}.vcf.gz".format(**locals())
		
		# Replace last colon (screws up file names)
		bcftools = """{xargs} --arg-file="{chrom_chunks_file}" -P {OPTIONS.cores} -I {{}} sh -c '{xarg_command}' """.format(**locals()).replace("__region__","{}")
		command(bcftools, c_log)
		
		#=========#
		# Filters #
		#=========#
		# Setup Filters here.

		# Heterozygous Polarization
		if COMMANDS.snps.bcftools.__heterozygous_polarization == True:
			filters = "| python {script_dir}/het_polarization.py ".format(**locals())
		else:
			filters = ""

		soft_filters = COMMANDS.snps.bcftools.__soft_filters
		hard_filters = COMMANDS.snps.bcftools.__hard_filters
		filters += construct_filters(COMMANDS.snps.bcftools.__soft_filters)
		print filters
		# Merge, sort, and index
		concat_list = ' '.join(["{vcf_dir}/TMP.{SM}.{chunk}.bcftools.{ind_union_filename}.vcf.gz".format(**locals()).replace(":","_") for chunk in chrom_chunks])
		bcftools_concat = """bcftools concat -O z {concat_list} {filters} > {complete_ind_vcf};
		 					 bcftools index {complete_ind_vcf}""".format(**locals()) 
		command(bcftools_concat, c_log)
		
		# Remove temporary files
		rm_comm = """rm {vcf_dir}/TMP.{SM}.*.bcftools.{ind_union_filename}.vcf.gz*""".format(**locals())
		command(rm_comm, c_log)



