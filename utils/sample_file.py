import csv
import yaml
import os
import glob
from pprint import pprint as pp
from utils import *

os.chdir("/Users/daniel/Documents/temp/")


class config:
    """
        Class for handling set up
    """
    def __init__(self, config):
        self.config = yaml.load(open(config, 'r'))
        script_dir = os.path.dirname(os.path.realpath(__file__)).replace("/utils", "")
        self.reqd_options = ["reference"
                             "fastq_dir",
                             "bam_dir",
                             "stat_dir",
                             "sample_file",
                             "chrom_chunk_kb"]

        # Set Options as base; for directories
        for k,i in self.config["OPTIONS"].items():
            setattr(self, k, i)
            if k.endswith("dir") and k != "fastq_dir":
                makedir(i)

        # Setup Reference Here
        ref_dir = "{script_dir}/genomes/{self.reference}/*fa.gz".format(**locals())
        self.reference_file = glob.glob(ref_dir)[0]


    def get_sample_file(self):
        return sample_file(self.sample_file, self)

    def chunk_genome(self, chunk_size):
        """ 
        Parses bwa .ann file to retrieve chromosome sizes
        for chunking purposes
        """
        ann = open(self.reference_file + ".ann").read()
        # Parsing .ann files
        contigs = [x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][::2]
        contig_sizes = map(int,[x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][1::2])
        chunk_size *= 1000
        for chrom, size in zip(contigs, contig_sizes):
            for chunk in xrange(1,size, chunk_size):
                if chunk + chunk_size > size:
                    chunk_end = size
                else:
                    chunk_end = chunk + chunk_size-1
                yield "{chrom}:{chunk}-{chunk_end}".format(**locals())


def get_non_uniq(non_uniq_list):
    non_uniq_list = list(set([x for x in non_uniq_list if non_uniq_list.count(x) > 1]))
    if len(non_uniq_list) > 0:
        return non_uniq_list
    else:
        return None

class sample_file:
    """
        Class for handling actions associated with the sample file.
    """

    def __init__(self, filename, config):
        self.sample_file = open(filename, 'rU')

        # Define Sets
        self.fastq_set = []
        self.bam_set = []
        self.vcf_set = []
        self.ID_set = []    # Used with Individual Bam Set.
        self.SM_ID_set = {} # Tracks bams (by ID) belonging to a sample.
        self.row_set = []

        self.required_values = ["FQ1", "FQ2", "ID", "LB", "SM"]
        with self.sample_file as f:
            iter_csv = csv.DictReader(f, delimiter = '\t', quoting = csv.QUOTE_NONE)
            for line, row in enumerate(iter_csv):
                # 1-based index
                line += 1 
                row["line"] = line

                # Track IDs
                ID = row["ID"]
                self.ID_set.append(ID)

                # set fq names.
                fq1, fq2 = row["FQ1"], row["FQ2"]
                fq_pair = tuple(config.fastq_dir + "/" + x for x in [fq1, fq2])
                # save fq pair
                row["fq_pair"] = fq_pair
                fq_exists = map(file_exists, fq_pair)
                self.fastq_set.append(fq_pair)

                empty_vals = [x for x in self.required_values if not row[x]]
                # Run basic checks
                if row["FQ1"] == row["FQ2"]:
                    raise Exception(
                        "Both Fastq's share same name on line %s: %s" % (line, fq1))
                elif row["SM"] == row["ID"]:
                    raise Exception(
                        "Sample Name and Sample ID can not be the same [line %s - %s]" % (line, row["ID"]))
                elif len(empty_vals) > 0:
                    empty_vals = ','.join(empty_vals)
                    raise Exception("Missing values on line %s: %s" % (line, empty_vals))
                elif not all(fq_exists):
                    missing_fastq = ','.join([fq_pair[x] for x in range(0,2) if not fq_exists[x]])
                    raise Exception(
                        "Fastq(s) Missing on line %s: %s" % (line, missing_fastq))

                if not row["PL"]:
                    row["PL"] = "ILLUMINA"

                # Construct a dictionary for read group for comparisons
                RG = dict([(k,v) for k,v in row.items() if k in ("ID", "LB", "SM", "PL",)])
                
                # Also Construct Raw Read Group (for aligning)
                rg_dat = [k + ":" + v for k,v in RG.items()]
                raw_RG = r"@RG\t" + r'\t'.join(sorted(rg_dat))

                row["RG"] = RG
                row["raw_RG"] = raw_RG

                # Group IDs and Samples by Read Group 
                SM = row["SM"]
                if SM not in self.SM_ID_set:
                    self.SM_ID_set[SM] = {"ID": [], "RG": []}
                self.SM_ID_set[SM]["ID"].append(ID)
                self.SM_ID_set[SM]["RG"].append(RG)

                # Remove keys incorporated into RG
                del row["LB"]
                del row["PL"]
                del row["SM"]

                # If all checks passed, append row dictionary
                self.row_set.append(row)

            # Sort Read Group List
            self.SM_ID_set
            # Check that all IDs are uniq
            non_uniq_IDs = get_non_uniq(self.ID_set)
            if get_non_uniq(self.ID_set):
                raise Exception(
                        "Non-uniq IDs exist: %s" % ', '.join(non_uniq_IDs))

            # Check that all fastq sets are uniq
            non_uniq_fastq_set = get_non_uniq(self.fastq_set)
            if non_uniq_fastq_set:
                raise Exception(
                        "Non-uniq Fastqs exist: %s" % ', '.join(non_uniq_fastq_set[0]))




config = config("/Users/daniel/Documents/temp/analysis.yaml")

sf = config.get_sample_file()
print pp(sf.SM_ID_set)
print pp(sf.row_set)
print pp(sf.SM_ID_set)
print [x for x in config.chunk_genome(100)]


