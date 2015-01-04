#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384
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

truncate_fq = """{stream_fq} {fq_loc} | head -n {OPTIONS.debug_fq_number_of_sequences} | gzip > {OPTIONS.fastq_dir}/DEBUG_FQ/{fq_filename}"""

bwa = """bwa mem -t {OPTIONS.cores} -R '{RG_header}' {bwa_options} {reference} {OPTIONS.fastq_dir}/{FQ1} {OPTIONS.fastq_dir}/{FQ2} | samtools view -@ {OPTIONS.cores} -bhu - > {bam_dir}/{ID}.unsorted.bam
         samtools -@ {OPTIONS.cores} sort -O bam -T {bam_dir}/{tmpname} {bam_dir}/{ID}.unsorted.bam > {bam_dir}/{ID}.sorted.bam"""

mark_dups = """
            java -jar {script_dir}/tools/picard.jar MarkDuplicates \
            I={bam_dir}/{ID}.sorted.bam \
            O={bam_dir}/{ID}.bam \
            M={bam_dir}/{ID}.duplicate_report.txt \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
            """

#====================#
# Load Configuration #
#====================#

opts = literal_eval(sys.argv[2])
config, log, c_log = load_config_and_log(sys.argv[1], "align")
OPTIONS = config.OPTIONS
COMMANDS = config.COMMANDS
align = COMMANDS.align # Pulls out alignment types.
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*fa.gz".format(**locals()))[0]
bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())
eav_file = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}/eav.txt".format(**locals())
stat_dir = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}".format(**locals())
eav = EAV()

#=============================#
# Setup Debug Mode if desired #
#=============================#
if OPTIONS.debug == True:
    # Truncate FASTQs for fast processing.
    for fq in ["FQ1", "FQ2"]:
        if not file_exists(OPTIONS.fastq_dir + "/DEBUG_FQ/" + opts[fq]):
            fq_loc = opts[fq.lower()]
            fq_filename = opts[fq]
            print fq_loc, fq_filename
            print truncate_fq.format(**locals())
            command(truncate_fq.format(**locals()), c_log)
    OPTIONS.fastq_dir = "{OPTIONS.fastq_dir}/DEBUG_FQ".format(**locals())
else:
    DEBUG = ""



#=========================#
# Setup Read Group Header #
#=========================#

# Set up Read Group String for alignment (with bwa)

ID = opts["ID"]
RG_header = construct_RG_header(ID, opts)


#===========================#
# Save Fastq Info and Stats #
#===========================#

fastq_info_file_loc = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}/FASTQ_INFO.txt".format(**locals())

cksum_set = []
try:
    cksum_set = Popen("grep 'cksum' {eav_file} | cut -f 5".format(**locals()), stdout=PIPE, shell=True).communicate()[0].strip().split("\n")
except:
    pass

opts["fq1_cksum"] = cksum(opts["fq1"])
opts["fq2_cksum"] = cksum(opts["fq2"])

# EAV remain throughout.
eav.file = eav_file
eav.entity = opts["SM"]


for fq in ["FQ1", "FQ2"]:
    # Check if fastq file stats already produced.
    if opts[fq.lower() + "_cksum"] not in cksum_set:
        eav.sub_entity = opts[fq]
        eav.attribute = "FASTQ Statistics".format(ID=opts["ID"])
        for k,v in extract_fastq_info(opts[fq.lower()]).items() + get_fastq_stats(opts[fq.lower()]).items():
            eav.sub_attribute = k
            eav.value = v
            eav.save()
        # Save Checksum
        eav.sub_attribute = "cksum"
        eav.value = opts[fq.lower() + "_cksum"]
        eav.save()
        # Save Read Group
        eav.sub_attribute = "Read Group"
        eav.value = RG_header


#=====#
# BWA #
#=====#

if "bwa" in align:
    bwa_options = format_command(align["bwa"])
    tmpname = os.path.split(tempfile.mktemp(prefix=ID))[1]
    FQ1 = opts["FQ1"]
    FQ2 = opts["FQ2"]

    # Create Directories
    makedir(bam_dir)
    completed_bam = "{bam_dir}/{ID}.bam".format(**locals())
    unsorted_bam = "{bam_dir}/{ID}.unsorted.bam".format(**locals())
    if not file_exists(completed_bam) and not file_exists(unsorted_bam):
        comm = bwa.format(**locals())
        command(comm, c_log)
        if align.alignment_options.remove_temp == True:
            file_to_delete = "rm {unsorted_bam}".format(**locals())
            command(file_to_delete, c_log)
    else:
        log.info("SKIPPING: " + completed_bam + " exists; no alignment.")


#=================#
# Mark Duplicates #
#=================#

def picard_dup_parser(dup_file, bam):
    histogram_values = ""
    with open(dup_file, 'r+') as f:
        for line in f:
            # Save aggregate stats
            if line.startswith("LIBRARY"):
                aggregate_values = bam + "\t" + f.next()
            # Save histogram stats
            elif line.startswith("BIN"):
                for i in range(0,100,1):
                    histogram_values += bam + "\t" + f.next()
    return (aggregate_values, histogram_values)


if "picard" in align:
    if "markduplicates" in align.picard:
        if align.picard.markduplicates == True:
            dup_report = "{bam_dir}/{ID}.duplicate_report.txt".format(**locals())
            if not file_exists(dup_report) or not file_exists(completed_bam):
                comm = mark_dups.format(**locals())
                log.info("Removing Duplicates: %s.bam" % ID)
                c_log.add(comm)
                command(comm, c_log)
                # Remove Sort tempfile
                if align.alignment_options.remove_temp == True:
                    file_to_delete = "rm {bam_dir}/{ID}.sorted.bam".format(**locals())
                    command(file_to_delete, c_log)
                # Process Duplicate Report
                stat_file = "{stat_dir}/PICARD_Aggregate_{OPTIONS.bam_dir}.txt".format(**locals())      # Leave OPTIONS.bam_dir
                histogram_file = "{stat_dir}/PICARD_Histogram_{OPTIONS.bam_dir}.txt".format(**locals()) # Leave OPTIONS.bam_dir
                # Aggregate Statistics from Picard dedup
                aggregate_values, histogram_values = picard_dup_parser(dup_report, ID + ".bam")
                if not file_exists(stat_file):
                    picard_stat_file = open(stat_file, 'w+')
                    picard_header = [ 'BAM',
                                      'LIBRARY',
                                      'UNPAIRED_READS_EXAMINED',
                                      'READ_PAIRS_EXAMINED',
                                      'UNMAPPED_READS',
                                      'ESTIMATED_LIBRARY_SIZE',
                                      'READ_PAIR_OPTICAL_DUPLICATES',
                                      'READ_PAIRS_EXAMINED', 'UNPAIRED_READS_EXAMINED',
                                      'PERCENT_DUPLICATION\n']
                    picard_stat_file.write("\t".join(picard_header))
                    picard_stat_file.close()
                with open(stat_file, "a") as stat:
                    stat.write(aggregate_values)
                if not file_exists(histogram_file):
                    picard_histogram_file = open(histogram_file,'w')
                    picard_histogram_file.write("BAM\tBIN\tVALUE\n")
                    picard_histogram_file.close()
                with open(histogram_file, 'a') as hist:
                    hist.write(histogram_values)
            else:
                log.info("SKIPPING: " + dup_report + " exists; Skipping.")
else:
    # If duplicates are not being marked, move files
    move_file = """mv {bam_dir}/{ID}.sorted.bam {bam_dir}/{ID}.bam""".format(**locals())
    command(move_file, c_log)

#=============================#
# Save BAM Stats (Individual) #
#=============================#

bam_individual = "{bam_dir}/{ID}.bam".format(**locals())
eav.sub_entity = ID + ".bam"
bam_individual_cksum = cksum(bam_individual)

if bam_individual_cksum not in cksum_set:
    for k,v in samtools_stats(bam_individual).items():
        print k,v
        eav.attribute = "BAM Statistics - Individual"
        eav.sub_attribute = k
        eav.value = v
        eav.save()
    eav.sub_attribute = "cksum"
    eav.value = bam_individual_cksum
    eav.save()




# Test for problems here..
sys.exit(0)
