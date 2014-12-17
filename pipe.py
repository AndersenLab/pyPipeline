#!/usr/bin/env python
"""pyPipeline

Usage:
  pipe.py align <config> [--debug]
  pipe.py samplefile <filename/dir>
  pipe.py genome [<name>]
  pipe.py test <config>

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
import glob
from utils import *
from utils.genomes import *
import csv

def check_fqs(fq):
    if not file_exists(fq["FQ1"]) or not file_exists(fq["FQ2"]):
        raise Exception("File Missing; Check: {fq1}, {fq2}".format(fq1=fq["FQ1"], fq2=fq["FQ2"]))

if __name__ == '__main__':
    opts = docopt(__doc__, version='pyPipeline')
    print opts

    #=======#
    # Setup #
    #=======#

    config_file = opts["<config>"]
    analysis_dir = os.getcwd()

    #==================#
    # Genome Retrieval #
    #==================#
    if opts["genome"] == True:
        if opts["<name>"] is not None:
            fetch_genome(opts["<name>"])
        else:
            list_genomes()
        exit()

    #=================#
    # New Sample File #
    #=================#

    """
        Create a sample file where Library, Sample, and Platform info can be added.
        Optionally update an analysis file
    """

    if opts["samplefile"] == True: 
        header = "FQ1\tFQ2\tLB\tSM\tPL\n"
        sample_file = open(opts["<filename/dir>"] + ".txt",'w')
        sample_file.write(header)
        if is_dir(analysis_dir + "/" + opts["<filename/dir>"]):
            # Construct a sample file using the directory info.
            fq_set = glob.glob(opts["<filename/dir>"] + "/*.fq.gz")
            fastq_pairs = zip(sorted([os.path.split(x)[1] for x in fq_set if x.find("1.fq.gz") != -1]), \
                sorted([os.path.split(x)[1] for x in fq_set if x.find("2.fq.gz") != -1]))
            for pair in fastq_pairs:
                sample_file.write("\t".join(pair) + "\n")
        exit()

    #===========#
    # Alignment #
    #===========#

    analysis_types = ["align", "merge", "snps","indels", "test"]
    analysis_type = [x for x in opts if opts[x] == True][0]
    # Load Configuration
    config, log, c_log = load_config_and_log(config_file, analysis_type)
    OPTIONS = config.OPTIONS
    log.info("#=== Beginning Analysis ===#")
    log.info("Running " + opts["<config>"])
    # Running locally or on a cluster
    if opts["--debug"] == True:
        run = "python"
        log.info("Using DEBUG mode")
    else:
        run = "sbatch"

    if analysis_type == "align":
        fq_set = open(config["OPTIONS"]["fastq_set"], 'rU')
        log.info("Performing Alignment")
        sample_set = {} # Generate a list of samples.
        for fq in csv.DictReader(fq_set, delimiter='\t', quoting=csv.QUOTE_NONE):
            fq1, fq2 = fq["FQ1"], fq["FQ2"]
            fq["FQ1"] = "{analysis_dir}/{OPTIONS.fastq_dir}/{fq1}".format(**locals())
            fq["FQ2"] = "{analysis_dir}/{OPTIONS.fastq_dir}/{fq2}".format(**locals())


            # Construct Individual BAM Dict
            SM = fq["SM"]
            if SM not in sample_set:
                sample_set[SM] = []
            ID = get_fq_ID([fq1, fq2])
            sample_set[SM].append(ID + ".bam")
            single_bam = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam".format(**locals())
            completed_merged_bam = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{SM}.bam".format(**locals())

            # Check that fq's exist before proceeding.
            check_fqs(fq)

            # Has the individually aligned bam been merged into a multi-bam bam?
            if file_exists(completed_merged_bam):
                pass
            # Does the individually aligned bam already exist?
            
            
            if not file_exists(single_bam):
                align = "{run} {script_dir}/align.py {config_file} \"{fq}\"".format(**locals())
                log.info(align)
                os.system(align)
            else:
                log.info("%-50s already aligned individually, skipping" % single_bam)
        # Merge Like Samples
        for SM in sample_set.keys():
            completed_merged_bam = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{SM}.bam".format(**locals())
            if not file_exists(completed_merged_bam):
                bams_to_merge = (SM, sample_set[SM])
                merge_bams = "{run} {script_dir}/merge_bams.py {config_file} \"{bams_to_merge}\"".format(**locals())
                log.info(merge_bams)
                os.system(merge_bams)
    if analysis_type == "test":
        reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*gz".format(**locals()))[0]
        for i in chunk_genome(3000000, reference):
            print i
        





    
