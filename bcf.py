# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, re, glob, hashlib, gzip, mimetypes, tempfile
from itertools import groupby as g
from collections import OrderedDict


def md5(filename):
    """Performs an md5 digest on a given file"""
    try:
        subprocess.check_output("which md5sum", shell=True) # Checks to see if md5sum exists.
        md5_command = "md5sum"
    except:
        md5_command = "md5"
    md5_result = subprocess.Popen('%s %s' % (md5_command, filename), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].replace("\n","")
    return md5_result.split("=")[1].strip()
 

def most_common(L):
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
 
class fastq(object):
    def __init__(self, fq1_filename, fq2_filename = None):
        self.fq1_filename = fq1_filename
        self.fq1 = self.extract_fastq_info(fq1_filename)
        self.fq_set = [self.fq1_filename]
        if fq2_filename is not None:
            self.fq2_filename = fq2_filename
            self.fq2 = self.extract_fastq_info(fq2_filename)
            self.fq_set += [self.fq2_filename]
        
    def extract_fastq_info(self, fq):
        """
        This function will extract information from the header lines of a demultiplexed fastq. Requires gzip to be imported.
        """
        line = {}
        if mimetypes.guess_type(fq)[1] == 'gzip':
            line["filetype"] = "gzip"
            f = gzip.open(fq, 'rb')
        else:
            line["filetype"] = "text"
            f = open(fq, 'rb')
        header_lines = [x.replace("\n","") for k,x in enumerate(f.readlines(10000)) if k%4==0]
        c = 0
        for heading in header_lines:
            c += 1
            l = re.split(r'(\:|#| )',heading)
            index_set = []
            if len(l) == 11:
                line["instrument"] = l[0]
                line["flowcell_lane"] = int(l[2])
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
                line["pair"] = int(l[14])
                line["filtered"] = l[16]
                line["control_bits"] = l[16]
                line["index"] = l[20]
                index_set.append(l[20])
            else:
                print "error", l
            line["index"] = most_common(index_set)
        line["fastq_filename"] = fq
        line["md5"] = md5(fq)
        return line
 
    def get_fastq_stats(self, fq_name):
        """
            This function will extract addition information from a fastq: Number of reads, unique reads, etc.
        """
        fq = getattr(self, fq_name)
        if fq["filetype"] == "gzip":
            command = "gunzip -c %s | awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'" % (fq["fastq_filename"], fq["fastq_filename"])
        else:
            command = "awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' %s" % (fq["fastq_filename"], fq["fastq_filename"])
        # Get stat results
        sr = subprocess.check_output(command, shell=True).strip().split(" ")
        fq["Number of Reads"] = sr[1]
        fq["Unique Reads"] = int(sr[2])
        fq["Frequency of Unique Reads"] = float(sr[3])
        fq["Most Abundant Sequence"] = sr[4]
        fq["Number of Times Most Abundant Sequence Occcurs"] = int(sr[5])
        fq["Frequency of Most Abundant Sequence"] = sr[6]

    def fastq_stats(self):
        self.get_fastq_stats("fq1")
        if len(self.fq_set) == 2:
            self.get_fastq_stats("fq2")

def format_options(options):
    formatted_options = {}
    # Force simple specification of some options for added functionality.
    simplified_ops = {"soft-filter": "-s"}
    for k,v in options.items():
        if k.replace("--","") in simplified_ops.keys():
            k = simplified_ops[k]
        if k.startswith("-"):
            formatted_options[k] = v
        elif len(k) == 1:
            formatted_options['-' + k] = v
        else:
            formatted_options['--' + k] = v
    # Wrap filters (include/exclude) in quotes
    for k,v in formatted_options.items():
        if k in ['--include','-i','--exclude','-e']:
            formatted_options[k] = "'%s'" % v
    return formatted_options

class bcf(file):
    def __init__(self, filename):
        # Start by storing basic information about the vcf/bcf and checking that an index exists.
        self.filename = filename
        self.filetype = os.path.splitext(filename)
        self.md5_digest = md5(filename)
        self.actions = []
        self.header_add_lines = []
        self.header = subprocess.Popen("bcftools view -h %s" % self.filename, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        # Fix index issues:
        old_index = self.header[1].startswith("Warning: The index file is older than the data file")
        no_index_file = os.path.isfile(filename + ".csi") == False
        if old_index or no_index_file and self.filetype in [".gz",".bcf"]:
            # Attempt to re-index if the index is older
            subprocess.check_output("bcftools index -f %s" % self.filename, shell=True)
            # Remove old indices.
            if os.path.isfile(self.filename + '.tbi'):
                os.remove(self.filename + ".tbi")
        self.header = self.header[0]

        # Samples
        self.samples = subprocess.check_output("bcftools query -l %s" % self.filename, shell=True).strip().split("\n")

        # Meta Data
        self.metadata = OrderedDict(re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M).findall(self.header))

        # Contigs
        self.contigs = [x.split(",") for x in re.compile(r'''^##contig=<(?P<data>.*)>''', re.M).findall(self.header)]
        self.contigs = [{x.split("=")[0]:x.split("=")[1] for x in f} for f in self.contigs]
        
        # Info
        self.info_set = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE).findall(self.header)
        
        # Filter
        self.filter_set = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE).findall(self.header)

        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE).findall(self.header)

    def variant_dict(self, lines=10):
        variant_set = []
        p = subprocess.Popen('bcftools view -H %s' % (self.filename), shell=True, stdout=subprocess.PIPE, bufsize=1)
        for line in iter(p.stdout.readline, b''):
            variant = {}
            line = line.strip().split("\t")
            variant["CHROM"] = line[0]
            variant["POS"] = int(line[1])
            variant["ID"] = line[2]
            variant["REF"] = line[3]
            variant["ALT"] = line[4]
            variant_set.append(variant)
            print variant

        p.communicate() 
        

    def filter(self, options, soft_filter_description=""):
        options["-O"] = "u" # Output in uncompressed bcf
        options = format_options(options)
        # If a soft-filter is used, add the description to the header.
        self.header_add_lines += ["##FILTER=<ID=%s,Description=\"%s\">" % (options["-s"], soft_filter_description)]
        options = ' '.join([k + " " + v for k,v in options.items()])
        self.actions += ["bcftools filter %s" % (options)]
        return self

    def stats(self):
        return subprocess.check_output("bcftools stats --samples - %s" % self.filename, stderr=subprocess.STDOUT, shell=True)

    def region(self, region):
        self.actions += ["bcftools view -O u -t {region}".format(region=region)]
        return self

<<<<<<< HEAD
=======
    def snpeff(self, annotation_db):
        # Apply snpeff annotations
        self.actions += ["bcftools view | snpeff eff %s" % annotation_db]
        self.header_add_lines += ["##INFO=<ID=EFF,Description=\"SNPEFF Annotation\">"]
        return self

>>>>>>> 713bf280620a5da0734ea6df4103391447cd7a89
    def out(self, out_filename, version=None):
        #===================#
        # Header Operations #
        #===================#

        self.header = self.header.split("\n")
        for l in self.header_add_lines:
            self.header.insert(1, l)

        # If version is 4.1, attempt to fix bcf/vcf so it is viewable in IGV.
        if version == 4.1:
            self.header[0] = "##fileformat=VCFv4.1"

        # Write out new header
        tmp_header = tempfile.NamedTemporaryFile()
        self.header = '\n'.join(self.header).strip()
        tmp_header.write(self.header)
        tmp_header.flush() # Flush out new header

        #====================================#
        # Add filters, reheader, annotations #
        #====================================#

        self.actions += ["bcftools reheader -h %s" % tmp_header.name]

<<<<<<< HEAD
        # Determine Output file type
        out_ext = os.path.splitext(out_filename)[1]
        out_opts = { ".bcf" : "b", ".gz" : "z", ".vcf" : "v"}
        if out_ext in out_opts.keys():
            self.actions += ["bcftools view -O %s" % out_opts[out_ext]]
        else:
            raise Exception("Unknown file extension (%s); Must Specify vcf, vcf.gz, or bcf" % out_filename)

=======
>>>>>>> 713bf280620a5da0734ea6df4103391447cd7a89
        if version == 4.1:
            self.actions += ["bcftools view | grep -v '##INFO' | bcftools view -O v"]

        # filetype
        file_output_types = {".bcf" : "b", ".gz" : "z", ".vcf" : "v"}
        filetype = file_output_types[os.path.splitext(out_filename)[1]]
        self.actions += ["bcftools view -O %s" % filetype]

        # Output types
        actions = ' | '.join(self.actions)

        print "bcftools view -O u %s | %s > %s" % (self.filename, actions, out_filename)
        subprocess.check_output("bcftools view -O u %s | %s > %s" % (self.filename, actions, out_filename), shell=True)
        if filetype in [".bcf", ".gz"]:
            subprocess.check_output("bcftools index %s" % out_filename, shell=True)

<<<<<<< HEAD
x = bcf("NIC276.nofilter.group.bcf")
x.filter({"include":'%QUAL>30', "soft-filter":"MaxQualityFail"})
x.filter({"include":'DP>3', "s": "MinimumDepth"})
x.out("NIC276.vcf.gz", version = 4.2)
=======
>>>>>>> 713bf280620a5da0734ea6df4103391447cd7a89

