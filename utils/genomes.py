import os, sys
import yaml
from utils import *

def fetch_genome(reference_name):
    """
        Downloads a reference genome and prepares for aligment
    """
    from utils import script_dir
    genome_list = yaml.load(open(script_dir + "/utils/genomes.yaml","r"))
    makedir("genomes")
    if reference_name not in genome_list:
        msg("Reference Genome not available", "error")
    ftp_loc = genome_list[reference_name]
    filename = os.path.split(ftp_loc)[1]
    makedir("{script_dir}/genomes/{reference_name}".format(**locals()))
    reference_loc = "{script_dir}/genomes/{reference_name}/{filename}".format(**locals())
    if not file_exists( reference_loc + ".sa"):
        print("Downloading {filename}".format(**locals()))
        os.system("curl {ftp_loc} > {script_dir}/genomes/{reference_name}/{filename}".format(**locals()))
        # Unzip and rezip with bgzip
        if filename.endswith(".gz"):
            os.system("gunzip {reference_loc} && bgzip {reference_loc2}".format(reference_loc=reference_loc, reference_loc2=reference_loc.replace(".gz","")))
        print("Indexing {script_dir}/genomes/{reference_name}/{filename}".format(**locals()))
        os.system("bwa index {script_dir}/genomes/{reference_name}/{filename}".format(**locals()))
    else:
        msg("Reference Already downloaded and indexed.", "error")

def list_genomes():
    """ 
        Prints a list of available genomes.
    """
    genome_list = yaml.load(open(script_dir + "/utils/genomes.yaml","r"))
    print("")
    print("\033[1m%-30s\t%-30s\033[0m" % ("Reference Name", "Location"))
    for k,v in genome_list.items():
        print("%-30s\t%-30s" % (k, v))
    print("")