#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
from commands import *
import tempfile
import glob

#====================#
# Load Configuration #
#====================#

opts = literal_eval(sys.argv[2])
config, log, c_log = load_config_and_log(sys.argv[1], "align")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
align = COMMANDS.align # Pulls out alignment types.

#=========================#
# Setup Read Group Header #
#=========================#

# Set up Read Group String for alignment (with bwa)
fqs = [os.path.split(opts["FQ1"])[1], os.path.split(opts["FQ2"])[1]]
ID = common_prefix(fqs).strip("-_")
SM = opts["SM"]
if opts["LB"] != "":
	LB = "LB:" + opts["LB"]
if opts["PL"] == "":
	PL = "PL:ILLUMINA"
else:
	PL = opts["PL"]
# Note library is optional; hence it's not explicitely defined.
RG_header = "@RG\\tID:{ID}\\tSM:{SM}\\t{PL}\\t{LB}".format(**locals())

#=====#
# BWA #
#=====#

if "bwa" in align:
	print config
	bwa_command, bwa_options = format_command(align["bwa"])
	reference = glob.glob("genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]
	tmpname = os.path.split(tempfile.mktemp(prefix=ID))[1]
	FQ1 = opts["FQ1"]
	FQ2 = opts["FQ2"]

	# Create Directories
	makedir(OPTIONS["analysis_dir"])
	makedir(OPTIONS["analysis_dir"] + "/bam")
	comm = bwa.format(**locals())
	command(comm, c_log)

#=================#
# Mark Duplicates #
#=================#

if "picard" in align:
	if "markduplicates" in align.picard:
		if align.picard.markduplicates == True:
			comm = mark_dups.format(**locals())
			log.info("Removing Duplicates: %s.bam" % ID)
			c_log.add(comm)
			command(comm, c_log)

