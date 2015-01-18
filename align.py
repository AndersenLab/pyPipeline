#!/usr/bin/python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384
import sys
import os
from ast import literal_eval
from utils import *
from utils.configuration import *
from utils.seq_utils import *
import tempfile
from pprint import pprint as pp

#=========#
# Command #
#=========#

truncate_fq = """{stream_fq} {fq_loc} | head -n {cf.debug_fq_number_of_sequences} | gzip > {cf.fastq_dir}/DEBUG_FQ/{fq_filename}"""

bwa = """bwa mem -t {cf.cores} -R '{RG_header}' {bwa_options} {reference} {cf.fastq_dir}/{FQ1} {cf.fastq_dir}/{FQ2} | samtools view -@ {cf.cores} -bhu - > {bam_dir}/{ID}.unsorted.bam
         samtools sort -@ {cf.cores} -O bam -T {bam_dir}/{tmpname} {bam_dir}/{ID}.unsorted.bam > {bam_dir}/{ID}.sorted.bam"""

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

job = dotdictify(literal_eval(sys.argv[2]))
cf = config(sys.argv[1])
sf = cf.get_sample_file()

#=============================#
# Setup Debug Mode if desired #
#=============================#
"""
if cf.debug == True:
    # Truncate FASTQs for fast processing.
    for fq in ["FQ1", "FQ2"]:
        if not file_exists(cf.fastq_dir + "/DEBUG_FQ/" + job[fq]):
            fq_loc = job[fq.lower()]
            fq_filename = job[fq]
            print fq_loc, fq_filename
            print truncate_fq.format(**locals())
            command(truncate_fq.format(**locals()), c_log)
    cf.fastq_dir = "{cf.fastq_dir}/DEBUG_FQ".format(**locals())
else:
    DEBUG = ""
"""

#===========================#
# Save Fastq Info and Stats #
#===========================#

cksum_set = []
try:
    cksum_set = Popen("grep 'cksum' {cf.eav.file} | cut -f 5".format(**locals()), stdout=PIPE, shell=True).communicate()[0].strip().split("\n")
except:
    pass

job.fq1_cksum = cksum(job.fq1)
job.fq2_cksum = cksum(job.fq2)
job.FQ1 = os.path.split(job.fq1)[1]
job.FQ2 = os.path.split(job.fq2)[1]

# EAV remain throughout.
eav = cf.eav
eav.entity = job["SM"]


for fq in ["FQ1", "FQ2"]:
    # Check if fastq file stats already produced.
    if job[fq.lower() + "_cksum"] not in cksum_set:
        print pp(job), "JOB"
        eav.sub_entity = job[fq]
        eav.attribute = "FASTQ Statistics".format(ID=job["ID"])
        for k, v in extract_fastq_info(job[fq.lower()]).items() + get_fastq_stats(job[fq.lower()]).items():
            eav.sub_attribute = k
            eav.value = v
            eav.save()
        # Save Checksum
        eav.sub_attribute = "cksum"
        eav.value = job[fq.lower() + "_cksum"]
        eav.save()
        # Save Read Group
        eav.sub_attribute = "Read Group"
        eav.value = job["@RG"]


#=====#
# BWA #
#=====#

if "bwa" in align:
    bwa_options = format_command(align["bwa"])
    tmpname = os.path.split(tempfile.mktemp(prefix=ID))[1]
    fq1_filename = os.path.split(job["fq1"])
    fq2_filename = os.path.split(job["fq2"])

    # Create Directories
    makedir(bam_dir)
    completed_bam = "{cf.bam_dir}/{job.ID}.bam".format(**locals())
    unsorted_bam = "{cf.bam_dir}/{job.ID}.unsorted.bam".format(**locals())
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
                stat_file = "{stat_dir}/PICARD_Aggregate_{cf.bam_dir}.txt".format(**locals())      # Leave cf.bam_dir
                histogram_file = "{stat_dir}/PICARD_Histogram_{cf.bam_dir}.txt".format(**locals()) # Leave cf.bam_dir
                # Aggregate Statistics from Picard dedup - seems error prone currently...
                try:
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
                except:
                    pass
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
        eav.attribute = "BAM Statistics - Individual"
        eav.sub_attribute = k
        eav.value = v
        eav.save()
    eav.sub_attribute = "cksum"
    eav.value = bam_individual_cksum
    eav.save()




# Test for problems here..
sys.exit(0)
