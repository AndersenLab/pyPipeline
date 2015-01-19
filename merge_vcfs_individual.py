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


#====================#
# Load Configuration #
#====================#

cf = config(sys.argv[1])
sf = cf.get_sample_file()
eav = cf.eav

#=========#
# Command #
#=========#

# Merging requires that filters within individual vcfs have passed.
merge_vcfs = """bcftools merge --apply-filters PASS -O z {vcf_set} > {vcf_dir}/{merged_vcf_name};
                bcftools index -f {vcf_dir}/{merged_vcf_name}"""


#=====================#
# Merge SNP VCF Files #
#=====================#

for merge_type in ["individual", "union"]:
    # If a union variant file already exists, skip individual merging
    print dir(cf)
    print cf.snps
    for caller in snp_callers:
        union_variant_file = "{OPTIONS.analysis_dir}/{caller}.{OPTIONS.union_variants}.txt".format(**locals())
        if (not file_exists(union_variant_file) and merge_type == "individual") or merge_type == "union":
            merged_vcf_name = merge_type + "_merged." + caller + ".vcf.gz"
            vcf_set = glob.glob("{vcf_dir}/*.{caller}.{merge_type}.vcf.gz".format(**locals()))
            vcf_set = ' '.join([x for x in vcf_set if not os.path.split(x)[1].startswith("individual") and not os.path.split(x)[1].startswith("TMP") and not os.path.split(x)[1].startswith("union")])
            if len(vcf_set) > 0:
                if not file_exists(vcf_set):
                    comm = merge_vcfs.format(**locals())
                    command(comm, c_log)
            # Generate a list of union variant sites among merged vcf
            if merge_type == "individual":
                ind_union_variant_set = "{OPTIONS.analysis_dir}/{caller}.{OPTIONS.union_variants}.txt".format(**locals())
                comm = """bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' {vcf_dir}/{merged_vcf_name} > {ind_union_variant_set}""".format(**locals())
                command(comm, c_log)
        if merge_type == "individual" and COMMANDS.snps.snp_options.remove_temp == True:
            glob.glob("{vcf_dir}/*{caller}.individual.vcf.gz*".format(**locals()))
            for i in glob.glob("{vcf_dir}/*{caller}.individual.vcf.gz*".format(**locals())):
                comm = "rm %s" % i
                cf.command(comm)