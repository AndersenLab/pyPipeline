import os, sys
import yaml
from utils import *

def fetch_genome(reference_name):
    """
        Downloads a reference genome and prepares for aligment
    """
    genome_list = yaml.load(open("utils/genomes.yaml","r"))
    makedir("genomes")
    ftp_loc = genome_list[reference_name]
    filename = os.path.split(ftp_loc)[1]
    makedir("genomes/{reference_name}".format(reference_name=reference_name))
    reference_loc = "genomes/{reference_name}/{filename}".format(**locals())
    if not file_greater_than_0( reference_loc + ".sa"):
        print("Indexing {filename}".format(**locals()))
        os.system("curl {ftp_loc} > genomes/{reference_name}/{filename}".format(**locals()))
        print("Indexing genomes/{reference_name}/{filename}".format(**locals()))
        os.system("bwa index genomes/{reference_name}/{filename}".format(**locals()))
    else:
        print("")
        print("Reference Already downloaded and indexed.")
        print("")

def list_genomes():
    """ 
        Prints a list of available genomes.
    """
    print os.getcwd()
    genome_list = yaml.load(open("utils/genomes.yaml","r"))
    print("")
    print("%-30s\t%-30s" % ("Reference Name", "Location"))
    for k,v in genome_list.items():
        print("%-30s\t%-30s" % (k, v))
    print("")