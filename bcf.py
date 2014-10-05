# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, re, glob, hashlib, gzip, mimetypes
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
    for k,v in options.items():
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
        self.md5_digest = md5(filename)
        self.actions = []
        self.header = subprocess.Popen("bcftools view -h %s" % self.filename, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if self.header[1].startswith("Warning: The index file is older than the data file"):
            # Attempt to re-index if the index is older.
            subprocess.check_output("bcftools index -f %s" % self.filename, shell=True)
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

    def filter(self, options):
        options["-O"] = "u" # Output in uncompressed bcf
        options = ' '.join([k + " " + v for k,v in format_options(options).items()])
        self.actions += ["bcftools filter %s" % (options)]
        return self

    def stats(self):
        return subprocess.check_output("bcftools stats --samples - %s" % self.filename, stderr=subprocess.STDOUT, shell=True)

    def region(self,chrom,start,end):
        print "bcftools filter -H -r %s:%s-%s %s" % (chrom, start, end, self.file)
        self.actions += ["bcftools view  -r %s:%s-%s %s" % (chrom, start, end, self.file)]
        return self

    def out(self, out_filename, type="bcf", version="4.1"):
        # If version is 4.1, attempt to fix bcf/vcf so it is viewable in IGV.

        if type is None:
            if out_filename.endswith(".bcf"):
                output_options = ""
                {'.bcf':"-O b", '.vcf' : "-O v", '.vcf.gz' : "-O z"}
        else:
            self.actions += ["bcftools view -O b"]
        # Output types
        actions = ' | '.join(self.actions)
        
        out_command = "bcftools view -O u %s | %s > %s" % (self.filename, actions, out_filename)
        print out_command
        return subprocess.check_output("bcftools view -O u %s | %s > %s" % (self.filename, actions, out_filename), shell = True)





x = fastq("BGI2-RET2-test1-f0534-1.fq.gz", "BGI1-RET2-test1-6b0f1-1.fq")

print x.fq1
print x.fastq_stats()
print x.fq2

#x = bcf("04_mmp_strains.txt.vcf.gz")
print x.filter({"include":'%QUAL>30', "soft-filter":"MaxQualityFail"})
print x.filter({"include":'DP>3', "soft-filter": "Minimum-Depth"})
print x.actions

x.out("fixed", "bcf")

