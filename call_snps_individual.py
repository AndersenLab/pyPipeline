#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
import tempfile
import glob


#====================#
# Load Configuration #
#====================#

opts = literal_eval(sys.argv[2])
config, log, c_log = load_config_and_log(sys.argv[1], "align")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
snps = COMMANDS.snps # Pulls out snp types.
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]
makedir("{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals()))

# Setup Chromosome Chunks
chrom_chunks = '\n'.join([x for x in chunk_genome(OPTIONS.chrom_chunk_kb*1000, reference)])
chrom_chunks_file = "{OPTIONS.analysis_dir}/chrom_chunks.txt".format(**locals())
with open(chrom_chunks_file, "w+") as f:
	f.write(chrom_chunks)
	
#==========#
# BCFTools #
#==========#

if 'bcftools' in snps:
	bam = "test1.bam"
	samtools_options = format_command(snps["samtools"])
	bcftools_options = format_command(snps["bcftools"])
	samtools_mpileup = "samtools mpileup {samtools_options} -g -f {reference} -r __region__ {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{bam}".format(**locals())
	bcftools_call = "bcftools call {bcftools_options}".format(**locals())
	xarg_command = "{samtools_mpileup} | {bcftools_call} > {OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}/{bam}.__region__.bcftools.bcf".format(**locals())
	bcftools = """gxargs --arg-file="{chrom_chunks_file}" -t -P {OPTIONS.cores} -I {{}} sh -c '{xarg_command}' """.format(**locals()).replace("__region__","{}")
	print bcftools
