#/usr/bin/python
import sys, os
from ast import literal_eval
from utils import *
import tempfile
import glob

#=========#
# Command #
#=========#

bwa = """bwa mem -R '{RG_header}' {bwa_options} {reference} {OPTIONS.fastq_dir}/{FQ1} {OPTIONS.fastq_dir}/{FQ2} | samtools view -bhu - > {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.unsorted.bam
         samtools sort -O bam -T {tmpname} {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.unsorted.bam > {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam"""

mark_dups = """
            java -jar {script_dir}/tools/picard.jar MarkDuplicates \
            I={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam \
            O={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam \
            M={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.duplicate_report.txt \
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
reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]

#=========================#
# Setup Read Group Header #
#=========================#

# Set up Read Group String for alignment (with bwa)
fqs = [os.path.split(opts["FQ1"])[1], os.path.split(opts["FQ2"])[1]]
ID = opts["ID"]
RG_header = construct_RG_header(ID, opts)

#=====#
# BWA #
#=====#

if "bwa" in align:
    bwa_options = format_command(align["bwa"])
    tmpname = os.path.split(tempfile.mktemp(prefix=ID))[1]
    FQ1 = opts["FQ1"]
    FQ2 = opts["FQ2"]

    # Create Directories
    makedir(OPTIONS["analysis_dir"])
    makedir(OPTIONS["analysis_dir"] + "/" + OPTIONS.bam_dir)
    makedir(OPTIONS["analysis_dir"] + "/statistics")
    completed_bam = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam".format(**locals())
    unsorted_bam = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.unsorted.bam".format(**locals())
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
            dup_report = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.duplicate_report.txt".format(**locals())
            if not file_exists(dup_report) or not file_exists(completed_bam):
                comm = mark_dups.format(**locals())
                log.info("Removing Duplicates: %s.bam" % ID)
                c_log.add(comm)
                command(comm, c_log)
                # Remove Sort tempfile
                if align.alignment_options.remove_temp == True:
                    file_to_delete = "rm {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam".format(**locals())
                    command(file_to_delete, c_log)
                # Process Duplicate Report
                stat_file = "{OPTIONS.analysis_dir}/statistics/PICARD_Aggregate_{OPTIONS.bam_dir}.txt".format(**locals())
                histogram_file = "{OPTIONS.analysis_dir}/statistics/PICARD_Histogram_{OPTIONS.bam_dir}.txt".format(**locals())
                # Aggregate Statistics from Picard dedup
                aggregate_values, histogram_values = picard_dup_parser(dup_report, ID + ".bam")
                if not file_exists(stat_file):
                    picard_stat_file = open(stat_file, "w+")
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
                with open(stat_file, "a+") as stat:
                    stat.write(aggregate_values)
                if not file_exists(histogram_file):
                    picard_histogram_file = open(histogram_file,"w+")
                    picard_histogram_file.write("BAM\tBIN\tVALUE\n")
                    picard_histogram_file.close()
                with open(histogram_file, 'a+') as hist:
                    hist.write(histogram_values)


            else:
                log.info("SKIPPING: " + dup_report + " exists; Skipping.")
else:
    # If duplicates are not being marked, move files
    move_file = """mv {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam""".format(**locals())
    command(move_file, c_log)
