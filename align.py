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

bwa = """bwa mem -t {cf.cores} -R '{job.RG}' {cf.align.bwa_options} {cf.reference_file} {fq1_filename} {fq2_filename} | samtools view -@ {cf.cores} -bhu - > {cf.bam_dir}/{job.ID}.unsorted.bam
         samtools sort -@ {cf.cores} -O bam -T {cf.bam_dir}/{tmpname} {cf.bam_dir}/{job.ID}.unsorted.bam > {cf.bam_dir}/{job.ID}.sorted.bam"""

mark_dups = """
            java -jar {script_dir}/tools/picard.jar MarkDuplicates \
            I={cf.bam_dir}/{job.ID}.sorted.bam \
            O={cf.bam_dir}/{job.ID}.bam \
            M={cf.bam_dir}/{job.ID}.duplicate_report.txt \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
            """

#====================#
# Load Configuration #
#====================#

job = dotdictify(literal_eval(sys.argv[2]))
cf = config(sys.argv[1])

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

if "bwa" in cf.align:
    tmpname = os.path.split(tempfile.mktemp(prefix=job.ID))[1]
    fq1_filename = job["fq1"]
    fq2_filename = job["fq2"]
    print job
    # Create Directories
    makedir(cf.bam_dir)
    completed_bam = "{cf.bam_dir}/{job.ID}.bam".format(**locals())
    unsorted_bam = "{cf.bam_dir}/{job.ID}.unsorted.bam".format(**locals())
    if not file_exists(completed_bam) and not file_exists(unsorted_bam):
        comm = bwa.format(**locals())
        cf.command(comm)
        if cf.align.alignment_options.remove_temp == True:
            file_to_delete = "rm {unsorted_bam}".format(**locals())
            cf.command(file_to_delete)
    else:
        cf.log("SKIPPING: " + completed_bam + " exists; no alignment.")


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


if cf.align.picard.markduplicates is True:
    dup_report = "{cf.bam_dir}/{job.ID}.duplicate_report.txt".format(**locals())
    if not file_exists(dup_report) or not file_exists(completed_bam):
        comm = mark_dups.format(**locals())
        cf.log("Removing Duplicates: %s.bam" % job.ID)
        cf.command(comm)
        # Remove Sort tempfile
        if cf.align.alignment_options.remove_temp is True:
            file_to_delete = "rm {cf.bam_dir}/{job.ID}.sorted.bam".format(**locals())
            cf.command(file_to_delete)
        # Process Duplicate Report
        stat_file = "PICARD_Aggregate_{cf.bam_dir}.txt".format(**locals())      # Leave cf.bam_dir
        histogram_file = "PICARD_Histogram_{cf.bam_dir}.txt".format(**locals()) # Leave cf.bam_dir
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
        cf.log("SKIPPING: " + dup_report + " exists; Skipping.")
else:
    # If duplicates are not being marked, move files
    move_file = """mv {cf.bam_dir}/{job.ID}.sorted.bam {cf.bam_dir}/{job.ID}.bam""".format(**locals())
    cf.command(move_file)

#=============================#
# Save BAM Stats (Individual) #
#=============================#

bam_individual = "{cf.bam_dir}/{job.ID}.bam".format(**locals())

# Index
cf.command("samtools index {bam_individual}".format(**locals()))

eav.sub_entity = job.ID + ".bam"
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
