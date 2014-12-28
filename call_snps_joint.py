#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
import tempfile
import glob
import csv

# Sent to node, calls snps by chromosome.

#====================#
# Load Configuration #
#====================#

chrom_region = sys.argv[2]
chrom_region_filename = chrom_region.replace(":", "_")
config, log, c_log = load_config_and_log(sys.argv[1], "snps")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
snps = COMMANDS.snps # Pulls out snp options.
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]

vcf_dir = "{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals())
bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())

makedir("{vcf_dir}".format(**locals()))

# Construct Bam List
sample_set = open(config.OPTIONS.sample_file, 'rU')
bam_set = []
for row in csv.DictReader(sample_set, delimiter='\t', quoting=csv.QUOTE_NONE):
    SM = row["SM"]
    bam_file = "{bam_dir}/{SM}.bam".format(**locals())
    bam_set.append(bam_file)
    if not file_exists(bam_file) or file_exists(bam_file + ".csi"):
        raise Exception("Bam File or index does not exist: %s" % bam_file)

bam_set = ' '.join(bam_set)

#==========#
# BCFTools #
#==========#
if 'bcftools' in snps:
	samtools_options = format_command(snps["samtools"])
	bcftools_options = format_command(snps["bcftools"])
	samtools_mpileup = "samtools mpileup {samtools_options} -g -f {reference} -r {chrom_region} {bam_set}".format(**locals())
	
	bcftools_call = "bcftools call --skip-variants indels {bcftools_options} ".format(**locals())
	
	#=========#
	# Filters #
	#=========#

	# Heterozygous Polarization
	if COMMANDS.snps.bcftools.__heterozygous_polarization == True:
		filters = "| python {script_dir}/het_polarization.py ".format(**locals())
	else:
		filters = ""

	soft_filters = COMMANDS.snps.bcftools.__soft_filters
	hard_filters = COMMANDS.snps.bcftools.__hard_filters
	filters += construct_filters(COMMANDS.snps.bcftools.__soft_filters)


	output_file = "{vcf_dir}/TMP.joint.{chrom_region_filename}.bcftools.vcf.gz".format(**locals())
	if not file_exists(output_file) or not file_exists(output_file + ".csi"):
		comm = "{samtools_mpileup} | {bcftools_call} {filters} > {output_file} && bcftools index {output_file}".format(**locals())
		os.system(comm)
	




