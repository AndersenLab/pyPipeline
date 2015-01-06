import os,sys
import re
import yaml
from subprocess import Popen, PIPE
import logging
from datetime import datetime



#======================#
# Set System Specifics #
#======================#

if os.uname()[0] == "Darwin":
    LOCAL = True
    xargs = "gxargs"
    stream_fq = "gunzip -kfc"
else:
    LOCAL = False
    xargs = "xargs"
    stream_fq = "zcat"

class dotdictify(dict):
    """
        http://stackoverflow.com/questions/3031219
        Creates dictionary accessible via attributes or keys.
    """
    marker = object()
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError, 'expected dict'

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, dotdictify):
            value = dotdictify(value)
        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        found = self.get(key, dotdictify.marker)
        if found is dotdictify.marker:
            found = dotdictify()
            dict.__setitem__(self, key, found)
        return found

    __setattr__ = __setitem__
    __getattr__ = __getitem__


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
        return "\nEntity:{self.entity}\nEntity:{self.sub_entity}\nAttribute:{self.attribute}\nSub-Attribute:{self.sub_attribute}\nValue:{self.value}\ntimestamp:{self.timestamp}\n".format(**locals())

    def save(self):
        if self.file is None:
            raise Exception("No Log File Set")
        if not file_exists(self.file):
            write_header = True
        else:
            write_header = False
        with(open(self.file,"a")) as f:
            if write_header == True:
                f.write("entity\tsub_entity\tattribute\tsub_attribute\tvalue\tcomment\ttimestamp\n")
            line = '\t'.join([self.entity, self.sub_entity, self.attribute, self.sub_attribute, str(self.value), self.comment, str(self.timestamp)])
            f.write(line.format(**locals()) + "\n")


def chunk_genome(chunk_size, reference):
    """ 
    Parses bwa .ann file to retrieve chromosome sizes
    for chunking purposes
    """
    ann = open(reference + ".ann").read()
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


def construct_filters(filter_list, soft=True):
    """ Constructs set of piped filters """
    if len(filter_list) > 0:
        filter_command = []
        for k,v in filter_list.items():
            if soft == True:
                filter_command.append("bcftools filter -O u --soft-filter {k} --exclude '{v}'".format(**locals()))
            else:
                filter_command.append("bcftools filter -O u --exclude '{v}'".format(**locals()))
        return '| ' + ' | '.join(filter_command) + ' | bcftools view -O z '
    else:
        return ''

def setup_logger(config):
    # Set up Logger
    log = logging.getLogger("pyPipeline")
    analysis_dir = config.OPTIONS.analysis_dir
    fh = logging.FileHandler(analysis_dir + "/" + analysis_dir + ".log")
    # Setup Formatting
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Set formats
    fh.setFormatter(formatter)
    log.addHandler(fh)
    if config.DEBUG == True:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    return log

class command_log:
    def __init__(self, config, job_type):
        analysis_dir = config.OPTIONS.analysis_dir
        self.log = open(analysis_dir + "/" + job_type + ".commands.log",'a')
    def add(self, command):
        # Clean up whitespace.
        command = re.sub("[^\S\r\n]+"," ", command).replace("\n ","\n").strip() + "\n"
        self.log.write(command)


def split_read_group(RG):
    RG_set = []
    for val in RG:
        val = sorted([x.split(":") for x in val.split("\t")[1:]])
        RG_set.append(val)
    return sorted(RG_set)

def load_config_and_log(config, job_type = None):
    """ 
        Loads the configuration file
            - Inherits default options not specified
            - Inherits options only for the job type specified.
    """
    default = dotdictify(yaml.load(open(script_dir + "/default.config.yaml","r")))
    config = dotdictify(yaml.load(open(config, 'r')))
    # Create analysis directory.
    makedir(config.OPTIONS.analysis_dir)
    # Override Options
    for opt,val in default["OPTIONS"].items():
        if opt not in config["OPTIONS"].keys():
            config["OPTIONS"][opt] = val
    general_log = setup_logger(config)
    c_log = command_log(config, job_type)

    if config.OPTIONS.debug == True:
        config.OPTIONS.analysis_dir = "DEBUG_" + config.OPTIONS.analysis_dir 

    return config, general_log, c_log

def format_command(command_config):
    """
        Performs standard formatting of commands being run.
    """
    opts = ""
    if command_config is not None:
        for k,v in command_config.items():
            # Use '__' for custom options
            if k.startswith("__"):
                pass
            # Use '_' to designate flags.
            elif k.startswith("_"):
                opts += " %s " % v
            else:
                opts += "%s %s " % (k,v)
        return opts
    else:
        return ""

def file_exists(filename):
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False

def makedir(dirname):
  if not os.path.exists(dirname):
    os.mkdir(dirname)

def is_dir(dir):
    return os.path.isdir(dir)

def remove_file(file):
    try:
        os.remove(file)
    except:
        pass

def command(command, log):
    """ Run a command on system and log """
    program = command.split(" ")[0]
    program_version = version(program)
    log.add("\n# " + program_version )
    log.add(command.strip() + "\n")
    command = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    for line in command.stdout:
        sys.write.stdout(line)
    if command.stderr is not None:
        raise Exception(command.stderr)

def which(program):
    """ returns the path to an executable or None if it can't be found
     http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
     """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def version(program):
    """ Attempts to return the version of an executable or None if it can't be found. """
    prog = which(program)
    if program not in ["tabix", "bwa", "rm", "cp", "ls"]:
        try:
            version = Popen([prog,"--version"], stdout=PIPE).communicate()[0]
            version = version.strip().split("\n")[0]
            if version == None:
                version = ""
            return "%-50s\t%s" % (prog, version)
        except:
            return ""
    else:
        # Hand special cases
        return "%-50s\t%s" % (prog, "")

def common_prefix(strings):
    """ Find the longest string that is a prefix of all the strings.
     http://bitbucket.org/ned/cog/src/tip/cogapp/whiteutils.py
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix

def get_script_dir():
    return os.path.dirname(os.path.realpath(__file__)).replace("/utils", "")

def is_defined(val):
    if val == '' or val == None:
        return False
    else:
        return True

def get_fq_ID(fqs):
    """ Returns common prefix of fastq's; stripping out select characters """
    return common_prefix(fqs).strip("-_")

def construct_RG_header(ID, opts):  
    if is_defined(opts["SM"]):
        SM = "SM:" + opts["SM"] + "\\t"
    else:
        SM = ""
    if is_defined(opts["LB"]):
        LB = "LB:" + opts["LB"] + "\\t"
    else:
        LB = ""
    if is_defined(opts["PL"]):
        PL = "PL:" + opts["PL"]
    else:
        PL = "PL:ILLUMINA"
    # Note library is optional; hence it's not explicitely defined.
    RG_header = "@RG\\tID:{ID}\\t{LB}{SM}{PL}".format(**locals())
    return RG_header

def get_bam_RG(bam):
    """
        Fetches Read Groups from the Header of SAM/BAM files
    """
    out, err = Popen(["samtools","view","-H",bam], stdout=PIPE).communicate()
    RG = [x for x in out.split("\n") if x.startswith("@RG")]
    return RG

def boolify(s):
    """ http://stackoverflow.com/questions/7019283/automatically-type-cast-parameters-in-python """
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def set_type(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def rreplace(s, old, new, count):
    """ 
    Replaces last occurance of something
    stackoverflow.com/questions/2556108/
    """
    return (s[::-1].replace(old[::-1], new[::-1], count))[::-1]

def get_column(filename, col_num, delim = "\t"):
    """ Returns list of a column from a file """
    column = []
    if file_exists(filename):
        with open(filename,'r') as r:
            for row in r:
                column.append(row.split(delim)[col_num])
        return column
    else:
        return []

# Define Constants
script_dir = get_script_dir()
available_snp_callers = ["bcftools", "freebayes"]

