#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
from commands import *
import tempfile
import glob

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

#================#
# Samtools Merge #
#================#

bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())

if len(SM_Bams) > 1:
    merge_options = format_command(align["merge"])
    merged_bam_name = SM + ".bam"
    SM_Bams = " ".join([bam_dir + "/" + x for x in SM_Bams])
    merge_bams = """samtools merge -f {merge_options} {bam_dir}/{merged_bam_name} {SM_Bams}""".format(**locals())
    command(merge_bams, c_log)
else:
    # Move single_sequence bam to merged name.
    ID = SM_Bams[0]
    move_file = """cp {bam_dir}/{ID} {bam_dir}/{SM}.bam""".format(**locals())
    command(move_file, c_log)

#================#
# Samtools Index #
#================#

command("samtools index {bam_dir}/{SM}.bam".format(**locals()), c_log)

#==================================#
# Remove individually aligned Bams #
#==================================#
if OPTIONS.alignment_options.remove_temp == True:
    for bam in SM_Bams:
        command("rm {bam_dir}/{bam}.bam".format(**locals()), c_log)
