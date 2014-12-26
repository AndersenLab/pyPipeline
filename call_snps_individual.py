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
makedir("{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals()))

# Setup Chromosome Chunks
chrom_chunks = [x for x in chunk_genome(OPTIONS.chrom_chunk_kb*1000, reference)]
chrom_list = '\n'.join([x for x in chunk_genome(OPTIONS.chrom_chunk_kb*1000, reference)])
chrom_chunks_file = "{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}_chrom_chunks.txt".format(**locals())
with open(chrom_chunks_file, "w+") as f:
	f.write(chrom_list)
	
#==========#
# BCFTools #
#==========#

if 'bcftools' in snps:
	samtools_options = format_command(snps["samtools"])
	bcftools_options = format_command(snps["bcftools"])
	samtools_mpileup = "samtools mpileup {samtools_options} -g -f {reference} -r __region__ {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{bam}".format(**locals())
	bcftools_call = "bcftools call --skip-variants indels {bcftools_options}".format(**locals())
	xarg_command = "REGION='__region__'; {samtools_mpileup} | {bcftools_call} -O z > {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/TMP.{SM}.${{REGION/:/_}}.bcftools.vcf.gz && bcftools index {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/TMP.{SM}.${{REGION/:/_}}.bcftools.vcf.gz".format(**locals())
	# Replace last colon (screws up file names)
	bcftools = """{xargs} --arg-file="{chrom_chunks_file}" -P {OPTIONS.cores} -I {{}} sh -c '{xarg_command}' """.format(**locals()).replace("__region__","{}")
	command(bcftools, c_log)
	# Merge, sort, and index
	concat_list = ' '.join(["{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/TMP.{SM}.{chunk}.bcftools.vcf.gz".format(**locals()).replace(":","_") for chunk in chrom_chunks])
	bcftools_concat = """bcftools concat -O z {concat_list} > {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/{SM}.bcftools.vcf.gz;
	 					 bcftools index {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/{SM}.bcftools.vcf.gz""".format(**locals()) 
	command(bcftools_concat, c_log)
	# Remove temporary files
	rm_comm = """rm {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/TMP.{SM}.*.bcftools.vcf.gz*""".format(**locals())
	command(rm_comm, c_log)
