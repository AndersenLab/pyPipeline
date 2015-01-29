import csv
import yaml
import os
import glob
from pprint import pprint as pp
from utils import *
from constants import *
from logs import *
from seq_utils import *
from subprocess import Popen, PIPE


class config:
    """
        Class for handling config set up.
    """
    def __init__(self, config):
        self.config_filename = config
        self.config_name = config.replace(".yaml", "")
        try:
            self.config = yaml.load(open(config, 'r'))
        except:
            msg("Sample File not found", "error")
        script_dir = os.path.dirname(os.path.realpath(__file__)).replace("/utils", "")

        self.node_index = 0

        # Set up Logs
        self.general_log = general_log(self)
        self.command_log = command_log(self)

        # Options that are required within the config file.
        self.reqd_options = ["reference",
                             "fastq_dir",
                             "bam_dir",
                             "log_dir",
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
                msg("You must define OPTION(s) %s in %s" % (undefined_options, analysis_filename), "error")

        # Set Options as base; for directories
        for k, i in self.config["OPTIONS"].items():
            setattr(self, k, i)
            if k.endswith("dir") and k != "fastq_dir":
                makedir(i)


        # Set Commands as base directories
        for analysis, values in self.config["COMMANDS"].items():
            setattr(self, analysis, dotdictify(values))
            for command, params in values.items():
                #setattr(self, command, dotdictify(params))
                if command in tools:
                    # Generate command options
                    opts = self.format_command_options(params)
                    analysis_attr = getattr(self, analysis)
                    setattr(analysis_attr, command + "_options", opts)
                    setattr(self, analysis, analysis_attr)


        # Setup command info
        self.cmd = dotdictify(self.config["COMMANDS"])

        # Set up entity-attribute-value
        self.eav = EAV()
        self.eav.file = self.config_name + ".eav.txt"

        # snp callers
        self.snp_callers = [x for x in self.cmd.snps if x in available_snp_callers]

        # Setup union variant sets (for SNP and INDEL)
        self.union_variants = dotdictify()
        for caller in self.snp_callers:
            for variant_type in ["ALL","SNP", "INDEL"]:
                self.union_variants[caller][variant_type] = "{self.config_name}.{variant_type}.{caller}.union_variants.txt".format(**locals())

        # Setup Reference Here
        ref_dir = "{script_dir}/genomes/{self.reference}/*f*.gz".format(**locals())
        try:
            self.reference_file = glob.glob(ref_dir)[0]
            assert(file_exists(self.reference_file))
        except:
            msg("Reference file '%s' does not exist" % (self.reference), "error")

    def log(self, msg, analysis_type=""):
        """ Adds to the log file for a given analysis. """
        self.general_log.add(msg, analysis_type)

    def command(self, command):
        """ Runs a command in the shell and logs that it was run. """

        out = Popen(command, shell=True, stdout=PIPE, stderr=None)
        for line in out.stdout:
            print(line)
        if out.stderr is not None:
            raise Exception(out.stderr)

    def format_command_options(self, command_config):
        """
            Performs standard formatting of commands being run.
        """
        opts = ""
        if command_config is not None:
            for k, v in command_config.items():
                # Use '__' for custom options
                if k.startswith("__"):
                    pass
                # Use '_' to designate flags.
                elif k.startswith("_"):
                    opts += " %s " % v
                else:
                    opts += "%s %s " % (k, v)
            return opts
        else:
            return ""

    def get_node():
        self.node_index += 1
        return str(self.nodes[node_index % len(self.nodes)])

    def submit_job(self,
                   command,
                   analysis_type,
                   log_name,
                   dependencies=None,
                   dependency_type="afterok"):
        """ Submit a job to the cluster """
        # Insert dependencies, output_dir, and nodes
        if LOCAL is False:
            self.log(command, "sbatch")
            command = command.split(" ")
            # Output Dirs
            log_file = "{self.log_dir}/{analysis_type}.{log_name}.%N.%j".format(**locals())
            output_dirs = " --output={log_file}.txt  --error={log_file}.err ".format(**locals())
            # Node
            if hasattr(self, "nodes"):
                use_node = "--nodelist={node} ".format(node="node" + get_node())
                command.insert(1, use_node)
            # Dependencies
            if dependencies is not None:
                if len(dependencies) > 0:
                    dependencies = ':'.join(dependencies)
                    depends_on = " --dependency={dep_type}:".format(**locals())
                    depends_on += dependencies
                else:
                    depends_on = ""
                command.insert(1, depends_on)
            else:
                depends_on = ""
            print command
            command = ' '.join(command)
            jobid, err = Popen(command, stdout=PIPE, stderr=PIPE, shell=True).communicate()
            jobid = jobid.strip().split(" ")[-1]
            print("Submitted job [{jobid}:{analysis_type}:{log_name}]".format(**locals()))
            if dependencies is not None:
                print("Dependencies: {dependencies}".format(dependencies=', '.join(dependencies)))

            if jobid.isdigit() is False:
                raise Exception("Error submitting %s" % jobid)
                exit()
            else:
                return jobid
        else:
            self.log(command, "python")
            self.command(command)

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
                chunk_set.append("{chrom}:{chunk}-{chunk_end}".format(**locals()).strip())
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
        self.filename = filename
        self.config = config
        self.sample_file_vars = ["FQ1", "FQ2", "ID", "LB", "SM"]
        # If the sample file exists, don't attempt to open it.
        if not file_exists(filename):
            return None
        self.sample_file = open(filename, 'rU')

        # Define Sets
        self.fastq_set = []
        self.ID_set = []        # Used with Individual Bam Set.
        self.SM_Group_set = dotdictify()  # Tracks bams (by ID) belonging to a sample.
        self.fq_set = []

        with self.sample_file as f:
            iter_csv = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            for line, row in enumerate(iter_csv):
                # 1-based index
                line += 1
                row["line"] = line

                # Track IDs
                ID = row["ID"]
                SM = row["SM"]
                self.ID_set.append(ID)

                # set fq names.
                fq1, fq2 = row["FQ1"], row["FQ2"]
                fq_pair = tuple(config.fastq_dir + "/" + x for x in [fq1, fq2])
                # save fq pair
                row["fq_pair"] = fq_pair
                fq_exists = map(file_exists, fq_pair)
                self.fastq_set.append(fq_pair)

                empty_vals = [x for x in self.sample_file_vars if not row[x]]
                # Run basic checks
                if row["FQ1"] == row["FQ2"]:
                    msg(
                        "Both Fastq's share same name on line %s in %s: %s" % (line, self.filename, fq1), "error")
                elif row["SM"] == row["ID"]:
                    msg(
                        "Sample Name and Sample ID can not be the same in %s [line %s - %s]" %
                        (self.filename, line, row["ID"]), "error")
                elif len(empty_vals) > 0:
                    empty_vals = ','.join(empty_vals)
                    msg("Missing values on line %s of %s: %s" % (line, self.filename, empty_vals), "error")
                elif not all(fq_exists):
                    missing_fastq = ','.join([fq_pair[x] for x in range(0, 2) if not fq_exists[x]])
                    msg(
                        "Fastq(s) Missing on line %s in %s: %s" % (line, self.filename, missing_fastq), "error")

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
                if SM not in self.SM_Group_set:
                    self.SM_Group_set[SM] = {"ID": [],
                                             "RG": [],
                                             "fq": [],
                                             "raw_RG": [],
                                             "bam_ind_filename": []}
                self.SM_Group_set[SM]["ID"].append(ID)
                self.SM_Group_set[SM]["RG"].append(RG)
                self.SM_Group_set[SM]["SM"] = SM
                self.SM_Group_set[SM]["RG"] = sorted(self.SM_Group_set[SM]["RG"])
                self.SM_Group_set[SM]["fq"].append(fq_pair)
                self.SM_Group_set[SM]["raw_RG"].append(raw_RG)
                self.SM_Group_set[SM]["bam_ind_filename"].append(bam_ind_filename)
                self.SM_Group_set[SM]["bam_merged_filename"] = bam_merged_filename

                # Add vcf files
                self.SM_Group_set[SM]["vcf_files"] = {}
                for caller in config.snp_callers:
                    for call_type in ["individual", "union"]:
                        for variant_type in ["ALL","SNP", "INDEL", "SV", "CNV", "TRANSPOSON"]:
                            vcf_ind = "{config.vcf_dir}/{SM}.{variant_type}.{caller}.{call_type}.vcf.gz".format(**locals())
                            self.SM_Group_set[SM]["vcf_files"][caller + "_" + call_type][variant_type] = vcf_ind


                # Remove keys incorporated into RG
                del row["LB"]
                del row["PL"]
                del row["SM"]

                # If all checks passed, append row dictionary
                self.fq_set.append(row)

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

    def new_sample_file(self):
        """
            Generates a new sample file from a directory
        """
        header = '\t'.join(self.sample_file_vars + ["RUN\n"])
        sample_filename = self.filename
        if is_dir(sample_filename):
            msg("Sample file is a directory", "error")
        new_sample_file = open(sample_filename, 'w')
        new_sample_file.write(header)
        sample_set = sorted(glob.glob(self.config.fastq_dir + "/*.fq.gz"))
        fastq_pairs = zip(sorted([os.path.split(x)[1] for x in sample_set if x.find("1.fq.gz") != -1]),
                          sorted([os.path.split(x)[1] for x in sample_set if x.find("2.fq.gz") != -1]))
        for pair in fastq_pairs:
            ID = common_prefix(pair).strip("-_")
            new_sample_file.write("\t".join(pair) + "\t" + ID + "\n")
        msg("Sample File Created")
        exit(0)

    def check_bams(self):
        """
            Generates list of bam merged (multi-sample)
            files; whether or not they exist, and whether
            or not they have the correct read group
            which reflects
        """
        for bam in self.SM_Group_set.values():
            bam["bam_merged_exists_and_RG_correct"] = False
            bam["bam_ind_exists_and_RG_correct"] = []
            if check_seq_file(bam["bam_merged_filename"]):
                if bam["RG"] == bamfile(bam["bam_merged_filename"]).RG:
                    bam["bam_merged_exists_and_RG_correct"] = True
                bam["bam_ind_exists_and_RG_correct"] = []
            for ind_bam in bam["bam_ind_filename"]:
                if check_seq_file(ind_bam):
                    if bamfile(ind_bam).RG[0] in bam["RG"]:
                        bam["bam_ind_exists_and_RG_correct"].append(True)
                else:
                    bam["bam_ind_exists_and_RG_correct"].append(False)
            yield dotdictify(bam)
