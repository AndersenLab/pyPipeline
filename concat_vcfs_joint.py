#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
from commands import *
import tempfile
import glob


#====================#
# Load Configuration #
#====================#

config, log, c_log = load_config_and_log(sys.argv[1], "align")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
vcf_dir = "{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals())

#=========#
# Command #
#=========#

# Merging requires that filters within individual vcfs have passed.
merge_vcfs = """bcftools concat -O z {vcf_list_string} > {merged_vcf_name};
                bcftools index -f {merged_vcf_name}"""


#=====================#
# Merge SNP VCF Files #
#=====================#

snp_callers = COMMANDS.snps
snp_callers = [x for x in snp_callers if x in available_snp_callers]
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]
chunks = [x.replace(":","_") for x in chunk_genome(OPTIONS.chrom_chunk_kb,reference)]

for caller in snp_callers:
    merged_vcf_name = "{vcf_dir}/joint.{caller}.vcf.gz".format(**locals())
    if not file_exists(merged_vcf_name) or not file_exists(merged_vcf_name + ".csi"):

        # If Merging - remove the index file to ensure temp files only deleted if things don't work.
        remove_file(merged_vcf_name + ".csi")

        vcf_list = ["{vcf_dir}/TMP.joint.{x}.{caller}.vcf.gz".format(**locals()) for x in chunks]
        # Check that all chunks are present.
        for vcf in vcf_list:
            if not file_exists(vcf) or not file_exists(vcf + ".csi"):
                raise Exception("VCF Chunk File or index missing: %s" % vcf)
        vcf_list_string = ' '.join(vcf_list)
        # Run Merge
        comm = merge_vcfs.format(**locals())
        command(comm, c_log)

        # If merge was successful, delete temporary files
        if file_exists(merged_vcf_name) and file_exists(merged_vcf_name + ".csi"):
            for vcf in vcf_list:
                comm = "rm {vcf} && rm {vcf}.csi".format(**locals())
                command(comm, c_log)

