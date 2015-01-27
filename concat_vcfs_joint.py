#!/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
from commands import *
import tempfile
import glob


#====================#
# Load Configuration #
#====================#

cf = config(sys.argv[1])

#=========#
# Command #
#=========#

# Merging requires that filters within individual vcfs have passed.
merge_vcfs = """bcftools concat -O z {vcf_list_string} > {merged_vcf_name};
                bcftools index -f {merged_vcf_name}"""


#=====================#
# Merge SNP VCF Files #
#=====================#

for caller in cf.snp_callers:
    joint_vcf_name = "{vcf_dir}/{cf.analysis_name}.{caller}.joint.vcf.gz".format(**locals())
    vcf_list = ["{vcf_dir}/TMP.joint.{x}.{caller}.vcf.gz".format(**locals()) for x in cf.chunk_genome()]
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

