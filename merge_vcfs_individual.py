#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=8192
import sys
import os
from utils import *
from utils.configuration import *
from commands import *
import glob
from pprint import pprint as pp

#====================#
# Load Configuration #
#====================#

cf = config(sys.argv[1])
sf = cf.get_sample_file()
eav = cf.eav

print pp(sf.SM_Group_set)

#=========#
# Command #
#=========#

# Merging requires that filters within individual vcfs have passed.
merge_vcfs = """bcftools merge --apply-filters PASS -O z {vcf_set} > {merged_vcf_name};
                bcftools index -f {merged_vcf_name}"""


#=====================#
# Merge SNP VCF Files #
#=====================#


# If a union variant file does not exist, perform individual calling.
for caller in cf.snp_callers:
    union_variant_file = cf.union_variants[caller]
    vcf_set = []
    if not file_exists(union_variant_file):
        merged_vcf_name = "{cf.vcf_dir}/{cf.config_name}.{caller}.individual.vcf.gz".format(**locals())
        for SM, data in sf.SM_Group_set.items():
            vcf_file = data["vcf_files"][caller + "_individual"]
            print vcf_file
            assert(file_exists(vcf_file))
            vcf_set.append(vcf_file)
        vcf_set = ' '.join(vcf_set)
        comm = merge_vcfs.format(**locals())
        cf.command(comm)
        # Run Individual
        comm = r"""bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' {merged_vcf_name} > {union_variant_file}""".format(**locals())
        cf.command(comm)
        # Remove individual
        if file_exists(union_variant_file) and all(check_seq_file(merged_vcf_name)):
            vcf_set = ' '.join([x for x in vcf_set.split(" ")] + [x + ".csi" for x in vcf_set.split(" ")])
            comm = "rm " + vcf_set
            cf.command(comm)
    else:
        print "RUNNING UNION"
        # Run Union
        merged_vcf_name = "{cf.vcf_dir}/{cf.config_name}.{caller}.union.vcf.gz".format(**locals())
        for SM, data in sf.SM_Group_set.items():
            vcf_file = data["vcf_files"][caller + "_union"]
            assert(file_exists(vcf_file))
            vcf_set.append(vcf_file)
        vcf_set = ' '.join(vcf_set)
        comm = merge_vcfs.format(**locals())
        cf.command(comm)