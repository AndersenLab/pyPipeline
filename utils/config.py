import csv
import yaml
import os
import glob
from pprint import pprint as pp
from utils import *
from seq_utils import *


os.chdir("/Users/Dan/Documents/tmp/")


class EAV:
    """
    Very simple Entity-Attribute-Value Object
    """

    def __init__(self):
        self.entity = ""
        self.sub_entity = ""
        self.attribute = ""
        self.sub_attribute = ""
        self.value = ""
        self.timestamp = datetime.now()
        self.comment = ""
        self.file = None

    def __repr__(self):
        return "\nEntity:{self.entity}\n\
                Entity:{self.sub_entity}\n\
                Attribute:{self.attribute}\n\
                Sub-Attribute:{self.sub_attribute}\n\
                Value:{self.value}\n\
                timestamp:{self.timestamp}\n".format(**locals())

    def save(self):
        if self.file is None:
            raise Exception("No Log File Set")
        if not file_exists(self.file):
            write_header = True
        else:
            write_header = False
        with(open(self.file, "a")) as f:
            if write_header is True:
                f.write("entity\tsub_entity\tattribute\tsub_attribute\tvalue\tcomment\ttimestamp\n")
            line = '\t'.join([self.entity,
                              self.sub_entity,
                              self.attribute,
                              self.sub_attribute,
                              str(self.value),
                              self.comment,
                              str(self.timestamp)])
            f.write(line.format(**locals()) + "\n")


class command_log:
    """
    The Command log is used to log commands that have been run.
    """
    def __init__(self, config, job_type):
        analysis_dir = config.OPTIONS.analysis_dir
        self.log = open(analysis_dir + "/" + job_type + ".commands.log", 'a')

    def add(self, command):
        # Clean up whitespace.
        command = re.sub("[^\S\r\n]+", " ", command).replace("\n ", "\n").strip() + "\n"
        self.log.write(command)


class config:
    """
        Class for handling configuration set up.
    """
    def __init__(self, config):
        self.config_filename = config
        self.config = yaml.load(open(config, 'r'))
        script_dir = os.path.dirname(os.path.realpath(__file__)).replace("/utils", "")

        # Options that are required within the config file.
        self.reqd_options = ["reference",
                             "fastq_dir",
                             "bam_dir",
                             "stat_dir",
                             "sample_file",
                             "chrom_chunk_kb",
                             "cores"]

        # Check that required options are defined
        for option in self.reqd_options:
            undefined_options = []
            for option in self.reqd_options:
                if option not in self.config["OPTIONS"].keys():
                    undefined_options.append(option)
            if undefined_options:
                undefined_options = ', '.join(undefined_options)
                analysis_filename = os.path.split(self.config_filename)[1]
                raise Exception("You must define OPTION(s) %s in %s" % (undefined_options, analysis_filename))

        # Set Options as base; for directories
        for k, i in self.config["OPTIONS"].items():
            setattr(self, k, i)
            if k.endswith("dir") and k != "fastq_dir":
                makedir(i)

        # Setup Reference Here
        ref_dir = "{script_dir}/genomes/{self.reference}/*f*.gz".format(**locals())
        try:
            self.reference_file = glob.glob(ref_dir)[0]
            assert(file_exists(self.reference_file))
        except:
            error("File does not exist")

    def get_sample_file(self):
        """
        Returns the sample file object
        """
        return sample_file(self.sample_file, self)

    def chunk_genome(self):
        """
        Parses bwa .ann file to retrieve chromosome sizes
        for chunking purposes
        """
        ann = open(self.reference_file + ".ann").read()
        # Parsing .ann files
        contigs = [x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][::2]
        contig_sizes = map(int, [x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][1::2])
        chunk_size = self.chrom_chunk_kb * 1000
        chunk_set = []
        for chrom, size in zip(contigs, contig_sizes):
            for chunk in xrange(1, size, chunk_size):
                if chunk + chunk_size > size:
                    chunk_end = size
                else:
                    chunk_end = chunk + chunk_size-1
                chunk_set.append("{chrom}:{chunk}-{chunk_end}".format(**locals()))
        return chunk_set


def get_non_uniq(non_uniq_list):
    non_uniq_list = list(set([x for x in non_uniq_list if non_uniq_list.count(x) > 1]))
    if len(non_uniq_list) > 0:
        return non_uniq_list
    else:
        return None


class sample_file:
    """
        Class for handling actions associated with the sample file:
    """

    def __init__(self, filename, config):
        self.sample_file = open(filename, 'rU')

        # Define Sets
        self.fastq_set = []
        self.bam_set = []
        self.vcf_set = []
        self.ID_set = []        # Used with Individual Bam Set.
        self.SM_Group_set = {}  # Tracks bams (by ID) belonging to a sample.
        self.fq_set = []

        self.required_values = ["FQ1", "FQ2", "ID", "LB", "SM"]
        with self.sample_file as f:
            iter_csv = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
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
                    missing_fastq = ','.join([fq_pair[x] for x in range(0, 2) if not fq_exists[x]])
                    raise Exception(
                        "Fastq(s) Missing on line %s: %s" % (line, missing_fastq))

                if not row["PL"]:
                    row["PL"] = "ILLUMINA"

                # Construct a dictionary for read group for comparisons
                RG = dict([(k, v) for k, v in row.items() if k in ("ID", "LB", "SM", "PL",)])

                # Also Construct Raw Read Group (for aligning)
                rg_dat = [k + ":" + v for k, v in RG.items()]
                raw_RG = r"@RG\t" + r'\t'.join(sorted(rg_dat))
                bam_ind_filename = config.bam_dir + "/" + row["ID"] + ".bam"
                bam_merged_filename = config.bam_dir + "/" + row["SM"] + ".bam"

                row["RG"] = RG
                row["raw_RG"] = raw_RG
                row["bam_ind_filename"] = bam_ind_filename
                row["bam_merged_filename"] = bam_merged_filename

                # Group IDs and Samples by Read Group
                SM = row["SM"]
                if SM not in self.SM_Group_set:
                    self.SM_Group_set[SM] = {"ID": [],
                                             "RG": [],
                                             "fq": [],
                                             "raw_RG": [],
                                             "bam_ind_filename": []}
                self.SM_Group_set[SM]["ID"].append(ID)
                self.SM_Group_set[SM]["RG"].append(RG)
                self.SM_Group_set[SM]["fq"].append(fq_pair)
                self.SM_Group_set[SM]["raw_RG"].append(raw_RG)
                self.SM_Group_set[SM]["bam_ind_filename"].append(bam_ind_filename)
                self.SM_Group_set[SM]["bam_merged_filename"] = bam_merged_filename

                # Remove keys incorporated into RG
                del row["LB"]
                del row["PL"]
                del row["SM"]

                # If all checks passed, append row dictionary
                self.fq_set.append(row)
                # Sort Read Group List
            self.SM_Group_set
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

    def get_merged_bams(self):
        """
            Generates list of bam merged (multi-sample)
            files; whether or not they exist, and whether
            or not they have the correct read group
            which reflects
        """
        print self.SM_Group_set

    def get_bam_individual(self):
        """
            Generates list of bam individual files,
            whether or not they have the correct
            read group, whether they exist, and
            whether or not their indices exist.
        """
        for row in self.fq_set:
            RG_from_sf = row["RG"]  # Read Group as defined in sample file.
            bam_ind_filename = row["bam_ind_filename"]
            bam_exists, bam_index_exists = check_seq_file(bam_ind_filename)
            RG_from_bam = False  # Define as false initially.
            if bam_exists:
                # Read Group as defined within aligned bam.
                RG_from_bam = bamfile(bam_ind_filename).RG[0]
            RG_identical = (RG_from_bam == RG_from_sf)
            yield {"RG_identical": RG_identical,
                   "bam_exists": bam_exists,
                   "bam_index_exists": bam_index_exists,
                   "bam_ind_filename": bam_ind_filename}






cf = config("analysis.yaml")

print dir(cf)

print cf.cores
print cf.reference_file

sf = cf.get_sample_file()

print pp(sf.SM_Group_set)
print pp(sf.fq_set)


#sf = config.get_sample_file()

#print pp([x for x in sf.get_merged_bams()])

