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
    config, gLOG, cLOG = load_config_and_logs(config_file, analysis_type)
    gLOG.info("Running " + opts["<config>"])
    # Running locally or on a cluster
    if opts["--debug"] == True:
        run = "python"
        gLOG.info("Using DEBUG mode")
    else:
        run = "sbatch"

    if analysis_type == "align":
        fq_set = open(config["OPTIONS"]["fastq_set"], 'rU')
        for fq in csv.DictReader(fq_set, delimiter='\t', quoting=csv.QUOTE_NONE):
            align = "{run} align.py {config_file} \"{fq}\"".format(**locals())
            command(align, cLOG)





    
