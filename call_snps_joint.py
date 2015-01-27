#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4096
import sys, os
from ast import literal_eval
from utils import *
from utils.configuration import *
import tempfile
import glob

# Sent to node, calls snps by chromosome.

#====================#
# Load Configuration #
#====================#

cf = config(sys.argv[1])
sf = cf.get_sample_file()
eav = cf.eav

chrom_region = sys.argv[2]
chunk_sanitized = chrom_region.replace(":","_")


bam_set = []
for i in sf.check_bams():
	bam_set.append(i["bam_merged_filename"])
bam_set = ' '.join(bam_set)

#==========#
# BCFTools #
#==========#
if 'bcftools' in cf.snps:
	samtools_mpileup = "samtools mpileup {cf.snps.samtools_options} -g -f {cf.reference_file} -r {chrom_region} {bam_set}".format(**locals())
	bcftools_call = "bcftools call --skip-variants indels {cf.snps.bcftools_options} ".format(**locals())

	#=========#
	# Filters #
	#=========#

	# Heterozygous Polarization
	if cf.snps.bcftools.__heterozygous_polarization == True:
		filters = "| bcftools view -O v | python {script_dir}/het_polarization.py ".format(**locals())
	else:
		filters = ""

	soft_filters = cf.snps.bcftools.__soft_filters
	hard_filters = cf.snps.bcftools.__hard_filters
	filters += construct_filters(cf.snps.bcftools.__soft_filters)

	output_file = "{cf.vcf_dir}/TMP.{cf.config_name}.joint.{chunk_sanitized}.bcftools.vcf.gz".format(**locals())
	if not file_exists(output_file) or not file_exists(output_file + ".csi"):
		comm = "{samtools_mpileup} | {bcftools_call} {filters} > {output_file} && bcftools index {output_file}".format(**locals())
		print comm
		os.system(comm)
	




