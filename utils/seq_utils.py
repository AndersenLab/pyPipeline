#!/usr/bin/python
import gzip
import re
from subprocess import Popen, PIPE
from itertools import groupby as g
from utils import *
import os


class bamfile:
    """
        Basic utility functions
    """

    def __init__(self, filename):
        self.header = self.fetch_header(filename)
        try:
            self.RG = self.header["RG"]
        except:
            self.RG = None

    def fetch_header(self, filename):
        comm = ["samtools", "view", "-H", filename]
        out, err = Popen(comm, stdout=PIPE).communicate()
        if err == "":
            raise Exception("Error parsing bam file.")
        else:
            header = out.strip().split("\n")
            header_set = {}
            for tag, val in [x.strip('@').split("\t", 1) for x in header]:
                if tag not in header_set:
                    header_set[tag] = []
                val = dict([tuple(x.split(':', 1)) for x in val.split('\t')])
                header_set[tag] += [val]
            header_set[tag] = sorted(header_set[tag])
            return header_set


def cksum(filename):
    """ Produces a file cksum """
    hash, err = Popen(["cksum", filename], stdout=PIPE, stderr=PIPE).communicate()
    if err != '':
        raise Exception("Error hashing {filename}".format(**locals()))
    return hash.split(" ")[0]


def most_common(L):
  try:
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
  except:
    return ""

def extract_fastq_info(fastq):
    """
    This function will extract information from the header lines of a demultiplexed fastq. Requires gzip to be imported.
    """
    f = gzip.open(fastq, 'rb')
    header_lines = [x.replace("\n","") for x in f.readlines(10000) if x.startswith("@")]

    for heading in header_lines:
            l = re.split(r'(\:|#| )',heading)
            line = {}
            index_set = []
            if len(l) == 11:
                line["instrument"] = l[0]
                line["flowcell_lane"] = l[2]
                line["flowcell_tile"] = l[4]
                try:
                    line["pair"] = l[10].split("/")[1]
                    index_set.append(l[10].split("/")[0])
                except:
                    pass
            elif len(l) == 21:
                line["instrument"] = l[0]
                line["run_id"] = l[2]
                line["flowcell_id"] = l[4]
                line["flowcell_lane"] = l[6]
                line["flowcell_tile"] = l[8]
                line["pair"] = l[14]
                line["filtered"] = l[16]
                line["control_bits"] = l[16]
                line["index"] = l[20]
                index_set.append(l[20])
            else:
                print "error", l
    line["index"] = most_common(index_set)
    return line


def get_fastq_stats(fq_name):
    """
        This function will extract addition information from a fastq: Number of reads, unique reads, etc.
    """
    fq = {}
    if fq_name.endswith("gz"):
        command = "gunzip -c %s | awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'" % (fq_name, fq_name)
    else:
        command = "awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' %s" % (fq_name, fq_name)
    # Get stat results
    sr = Popen(command, stdout=PIPE, shell=True).communicate()[0].strip().split(" ")
    fq["Number_of_Reads"] = sr[1]
    fq["Unique_Reads"] = int(sr[2])
    fq["Frequency_of_Unique_Reads"] = float(sr[3])
    fq["Most_Abundant_Sequence"] = sr[4]
    fq["Number_of_Times_Most_Abundant_Sequence_Occurs"] = int(sr[5])
    fq["Frequency_of_Most_Abundant_Sequence"] = sr[6]
    return fq


def get_contigs(bam):
    header, err = Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    # Extract contigs from header and convert contigs to integers
    contigs = {}
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
        contigs[x[0]] = int(x[1])
    return contigs

def coverage(bam, mtchr = None):
    # Check to see if file exists
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")
    contigs = get_contigs(bam)

    # Guess mitochondrial chromosome
    mtchr = [x for x in contigs if x.lower().find("m") == 0]
    if len(mtchr) != 1:
        mtchr = None
    else:
        mtchr = mtchr[0]

    coverage_dict = {}
    for c in contigs.keys():
        command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
        coverage_dict[c] = {}
        try:
            coverage_dict[c]["Bases Mapped"], coverage_dict[c]["Sum of Depths"] = map(int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))
            coverage_dict[c]["Breadth of Coverage"] = coverage_dict[c]["Bases Mapped"] / float(contigs[c])
            coverage_dict[c]["Depth of Coverage"] = coverage_dict[c]["Sum of Depths"] / float(contigs[c])
            coverage_dict[c]["Length"] = int(contigs[c])
        except:
            coverage_dict[c]["Bases Mapped"] = 0
            coverage_dict[c]["Sum of Depths"] = 0
            coverage_dict[c]["Breadth of Coverage"] = 0
            coverage_dict[c]["Depth of Coverage"] = 0
            coverage_dict[c]["Length"] = 0


    # Calculate Genome Wide Breadth of Coverage and Depth of Coverage
    genome_length = float(sum(contigs.values()))
    coverage_dict["genome"] = {}
    coverage_dict["genome"]["Length"] = int(genome_length)
    coverage_dict["genome"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    coverage_dict["genome"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)

    if mtchr != None:
        # Calculate nuclear breadth of coverage and depth of coverage
        ignore_contigs = [mtchr, "genome", "nuclear"]
        coverage_dict["nuclear"] = {}
        coverage_dict["nuclear"]["Length"] = sum([x["Length"] for k,x in coverage_dict.iteritems() if k not in ignore_contigs ])
        coverage_dict["nuclear"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])
        coverage_dict["nuclear"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])

        # Calculate the ratio of mtDNA depth to nuclear depth
        coverage_dict["genome"]["mt_ratio"] = coverage_dict[mtchr]["Sum of Depths"] / float(coverage_dict["nuclear"]["Depth of Coverage"])

    # Flatten Dictionary 
    coverage = []
    for k,v in coverage_dict.items():
        for x in v.items():
            coverage += [(x[0],k, x[1])]
    return coverage

def samtools_stats(filename):
    """ Load samtools stats into python """
    stats, err = Popen(["samtools","stats",filename], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    stats = [x.split("\t") for x in stats.split("\n")]
    chksum = [x for x in stats if x[0].startswith("CHK")][0]
    stats = dict([(x[1].replace(":",""),set_type(x[2]),) for x in stats if x[0].startswith("SN")])
    stats["filename"] = filename
    stats["chksum_read_names"] = chksum[1]
    stats["chksum_sequences"] = chksum[2]
    stats["chksum_qualities"] = chksum[3]
    return stats


def check_seq_file(filename):
    """
        Check whether a BAM or VCF and
        its associated indices exists
    """
    if filename.endswith(".vcf.gz") or filename.endswith(".bcf"):
        return all([file_exists(filename), file_exists(filename + ".csi")])
    if filename.endswith(".bam"):
        return all([file_exists(filename), file_exists(filename + ".bai")])


