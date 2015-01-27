#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=8192
import sys, os
from ast import literal_eval
from utils import *
from utils.seq_utils import *
from utils.configuration import *
from commands import *
import tempfile
import glob
from pprint import pprint as pp
# ======= #
# Command #
# ======= #

merge_bams = """sambamba merge -t {cf.cores} {merge_options} {bam_dir}/{merged_bam_name} {bam_dir}/{SM_Bams}"""

#====================#
# Load Configuration #
#====================#

job = literal_eval(sys.argv[2])
SM = job["SM"]
SM_Bams = job["bam_ind_filename"]
cf = config(sys.argv[1])
eav = cf.eav

merged_bam_name = cf.bam_dir + "/" + SM + ".bam"


#================#
# Samtools Merge #
#================#

bam_dir = "{cf.bam_dir}".format(**locals())
move_file = False


if len(SM_Bams) > 1:
    SM_Bam_set = " ".join(SM_Bams)
    merge_bams = """samtools merge -f {merged_bam_name} {SM_Bam_set}""".format(**locals())
    print merge_bams
    cf.command(merge_bams)
    # Index
    cf.command("samtools index {merged_bam_name}".format(**locals()))
else:
    # Move single_sequence bam to merged name.
    ID = SM_Bams[0]
    move_file = """mv {ID} {merged_bam_name};
                   mv {ID}.bai {merged_bam_name}.bai""".format(**locals())
    cf.command(move_file)
    dup_report = ID.replace(".bam",".duplicate_report.txt")
    cf.command("rm {dup_report}".format(**locals()))

#========================#
# Remove temp files here #
#========================#
if cf.align.alignment_options.remove_temp is True and move_file is False:
    for bam_ind in SM_Bams:
        cf.command("rm {bam_ind} && rm {bam_ind}.bai".format(**locals()))
        dup_report = bam_ind.replace(".bam",".duplicate_report.txt")
        cf.command("rm {dup_report}".format(**locals()))

#=======================#
# Collect Stats for BAM #
#=======================#

cksum_set = []
try:
    cksum_set = Popen("grep 'cksum' {cf.eav_file} | grep 'BAM Statistics - Merged' | cut -f 5".format(**locals()), stdout=PIPE, shell=True).communicate()[0].strip().split("\n")
except:
    pass

eav.entity = SM
eav.sub_entity = SM + ".bam"
eav.attribute = "BAM Statistics - Merged"
bam_merged_cksum = cksum(merged_bam_name)

if bam_merged_cksum not in cksum_set:
    # Save coverage info
    for contig, k,v in coverage(merged_bam_name):
        eav.sub_attribute = contig + " (" + k + ")"
        eav.value = v
        eav.save()
    for k,v in samtools_stats(merged_bam_name).items():
        eav.sub_attribute = k
        eav.value = v
        eav.save()
    eav.sub_attribute = "cksum"
    eav.value = bam_merged_cksum 
    eav.save()


# Test for problems here...
sys.exit(0)
