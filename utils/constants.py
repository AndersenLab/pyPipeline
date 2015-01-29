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
         "samtools",
         "freebayes",
         "bcftools",
         "picard"]

if os.uname()[0] == "Darwin":
    LOCAL = True
    echo = "gecho"
    xargs = "gxargs"
    sort = "gsort"
    run = "python"
    output_dirs = ""
    stream_fq = "gunzip -kfc"
else:
    run = "sbatch"
    echo = "echo"
    sort = "sort"
    LOCAL = False
    xargs = "xargs"
    stream_fq = "zcat"
