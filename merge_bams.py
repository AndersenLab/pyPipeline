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
from commands import *
import tempfile
import glob
import operator

# ======= #
# Command #
# ======= #

merge_bams = """samtools merge {merge_options} {bam_dir}/{merged_bam_name} {bam_dir}/{SM_Bams}"""

#====================#
# Load Configuration #
#====================#

SM, SM_Bams = literal_eval(sys.argv[2])
config, log, c_log = load_config_and_log(sys.argv[1], "align")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
align = COMMANDS.align # Pulls out alignment types.

bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())
eav_file = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}/eav.txt".format(**locals())
stat_dir = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}".format(**locals())
eav = EAV()
eav.file = eav_file

#================#
# Samtools Merge #
#================#

bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())

if len(SM_Bams) > 1:
    merge_options = format_command(align["merge"])
    merged_bam_name = SM + ".bam"
    SM_Bam_set = " ".join([bam_dir + "/" + x for x in SM_Bams])
    merge_bams = """samtools merge -f {merge_options} {bam_dir}/{merged_bam_name} {SM_Bam_set}""".format(**locals())
    command(merge_bams, c_log)
    # Index
    command("samtools index {bam_dir}/{SM}.bam".format(**locals()), c_log)
else:
    # Move single_sequence bam to merged name.
    ID = SM_Bams[0]
    move_file = """mv {bam_dir}/{ID} {bam_dir}/{SM}.bam;
                   mv {bam_dir}/{ID}.bai {bam_dir}/{SM}.bam.bai""".format(**locals())
    command(move_file, c_log)




#========================#
# Remove temp files here #
#========================#

if COMMANDS.align.alignment_options.remove_temp == True:
    print SM_Bams
    for bam in SM_Bams:
        command("rm {bam_dir}/{bam}".format(**locals()), c_log)
        dup_report = bam.replace(".bam",".duplicate_report.txt")
        command("rm {bam_dir}/{dup_report}".format(**locals()), c_log)

#=======================#
# Collect Stats for BAM #
#=======================#

cksum_set = []
try:
    cksum_set = Popen("grep 'cksum' {eav_file} | grep 'BAM Statistics - Merged' | cut -f 5".format(**locals()), stdout=PIPE, shell=True).communicate()[0].strip().split("\n")
except:
    pass

bam_merged = "{bam_dir}/{SM}.bam".format(**locals())
eav.entity = SM
eav.sub_entity = SM + ".bam"
eav.attribute = "BAM Statistics - Merged"
bam_merged_cksum = cksum(bam_merged)

if bam_merged_cksum not in cksum_set:
    # Save coverage info
    for contig, k,v in coverage(bam_merged):
        eav.sub_attribute = contig + " (" + k + ")"
        eav.value = v
        eav.save()
    for k,v in samtools_stats(bam_merged).items():
        eav.sub_attribute = k
        eav.value = v
        eav.save()
    eav.sub_attribute = "cksum"
    eav.value = bam_merged_cksum 
    eav.save()


# Test for problems here...
sys.exit(0)
