#!/usr/bin/env python
"""pyPipeline

Usage:
  pipe.py align <config> [--debug]
  pipe.py new (samplefile) <filename>
  pipe.py genome [<name>]

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
from utils import *
from utils.genomes import *
import csv

if __name__ == '__main__':
    opts = docopt(__doc__, version='pyPipeline')
    print opts

    #==================#
    # Genome Retrieval #
    #==================#
    if opts["genome"] == True:
        if opts["<name>"] is not None:
            fetch_genome(opts["<name>"])
        else:
            list_genomes()
        exit()

    #===========#
    # Alignment #
    #===========#
    if opts["new"] == True: 
        sample_file = "FQ1\tFQ2\tLB\tSM\tPL"
        open(opts["<filename>"],'w').write(sample_file)
        exit()
    analysis_types = ["align", "snps","indels"]
    analysis_type = [x for x in opts if opts[x] == True][0]
    # Load Configuration
    config_file = opts["<config>"]
    config, log, c_log = load_config_and_log(config_file, analysis_type)
    print(config)
    log.info("#============== Beginning Analysis ==============#")
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
        for fq in csv.DictReader(fq_set, delimiter='\t', quoting=csv.QUOTE_NONE):
            align = "{run} {script_dir}align.py {config_file} \"{fq}\"".format(**locals())
            log.info(align)
            os.system(align)
        # Merge Like Samples
        





    
