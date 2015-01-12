import sys, os
from ast import literal_eval
from utils import *
from utils.seq_utils import *
import tempfile
import glob
import csv

#=========#
# Command #
#=========#

mcclintock = """
              {run} mcclintock.sh
              -r {cf.reference_file}" \
              -c {consensus} \
              -1 {cf.fastq_dir}/{FQ1}\
              -2 {cf.fastq_dir}/{FQ2}\
              """
 


#====================#
# Load Configuration #
#====================#


#job = dotdictify(literal_eval(sys.argv[2]))

#cf = config(sys.argv[1])
cf = "analysis.yaml"
job = {"fq1": }
sf = cf.get_sample_file()



#change to locals
#commands
#subprocess.Popen("sh mcclintock.sh -r %s-c %s-1 %s-2 %s", %(reference, consensus, fastq1, fastq2)
subprocess.Popen(mcclintock.format(**locals()))

comm = mcclintock.format(**locals())
command(comm, transposons.command.log)
