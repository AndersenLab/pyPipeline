#!/usr/bin/python
import sys
import os
from ast import literal_eval
from utils import *
from utils.configuration import *
from utils.seq_utils import *
import tempfile
from pprint import pprint as pp


#====================#
# Load Configuration #
#====================#

cf = config(sys.argv[1])
sf = cf.get_sample_file()

#=========#
# Command #
#=========#

# Merging requires that filters within individual vcfs have passed.
merge_vcfs = """bcftools concat -O z {vcf_list_string} > {joint_vcf_name};
                bcftools index -f {joint_vcf_name}"""

#=====================#
# Merge SNP VCF Files #
#=====================#

for caller in cf.snp_callers:
    joint_vcf_name = "{cf.vcf_dir}/{cf.config_name}.{caller}.joint.vcf.gz".format(**locals())
    vcf_list = []
    for chunk in cf.chunk_genome():
        # Check that chunk does not exist.
        chunk_sanitized = chunk.replace(":","_")
        vcf_list.append("{cf.vcf_dir}/TMP.{cf.config_name}.joint.{chunk_sanitized}.{caller}.vcf.gz".format(**locals()))
    # Check that all chunks are present.
    for vcf in vcf_list:
        if not file_exists(vcf) or not file_exists(vcf + ".csi"):
            raise Exception("VCF Chunk File or index missing: %s" % vcf)
    vcf_list_string = ' '.join(vcf_list)
    # Run Merge
    comm = merge_vcfs.format(**locals())
    cf.command(comm)

    # If merge was successful, delete temporary files
    if all(check_seq_file(joint_vcf_name)):
        for vcf in vcf_list:
            comm = "rm {vcf} && rm {vcf}.csi".format(**locals())
            cf.command(comm)

