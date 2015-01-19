# Define constants
import os

script_dir = os.path.dirname(os.path.realpath(__file__)).replace("/utils", "")
available_snp_callers = ["bcftools", "freebayes"]
analysis_types = ["trim",
                  "align",
                  "merge",
                  "snps",
                  "indels",
                  "transposons",
                  "test"]

tools = ["bwa",
         "freebayes",
         "bcftools",
         "picard"]

if os.uname()[0] == "Darwin":
    LOCAL = True
    xargs = "gxargs"
    run = "python"
    output_dirs = ""
    stream_fq = "gunzip -kfc"
else:
    run = "sbatch"
    LOCAL = False
    xargs = "xargs"
    stream_fq = "zcat"

