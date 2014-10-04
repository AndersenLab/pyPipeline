# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, re, glob, hashlib, gzip
from collections import OrderedDict


def md5(filename):
    """Performs an md5 digest on a given file"""
    try:
        return subprocess.check_output(['md5sum',filename], shell=True).replace("\n","")
    except:
        return subprocess.check_output(['md5',filename], shell=True).replace("\n","")
 

def most_common(L):
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
 
class fastq(object):
    def __init__(self, fq1_filename, fq2_filename = None):
        self.fq1_filename = fq1_filename
        self.fq1 = self.extract_fastq_info(fq1_filename)
        if fq2_filename is not None:
            self.fq2_filename = fq2_filename
            self.fq2 = self.extract_fastq_info(fq2_filename)
 
    def extract_fastq_info(self, fq):
        """
        This function will extract information from the header lines of a demultiplexed fastq. Requires gzip to be imported.
        """
        f = gzip.open(fq, 'rb')
        header_lines = [x.replace("\n","") for x in f.readlines(10000) if x.startswith("@")]
 
        for heading in header_lines:
                l = re.split(r'(\:|#| )',heading)
                line = {}
                index_set = []
                if len(l) == 11:
                    line["instrument"] = l[0]
                    line["flowcell_lane"] = l[2]
                    line["flowcell_tile"] = l[4]
                    line["x_coord"] = l[6]
                    line["y_coord"] = l[8]
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
                    line["x_coord"] = l[10]
                    line["y_coord"] = l[12]
                    line["pair"] = l[14]
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
 
    def fastq_stats(self):
        fq_set = [self.fq1_filename, self.fq2_filename]
        try:
            # If the user has gnu parallel, both fastqs can be run at the same time.
            command = "parallel \"gunzip -c {} | awk -v filename={} '((NR-2)%%4==0){read=\$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'\" ::: %s" % ' '.join(fq_set)
            # Parse Stats
            fqstats = [x.split("\t") for x in subprocess.check_output(command, shell=True).strip().split("\n")]
            self.fq_statistics = {
                "Number of Reads": 1
            }
        except:
            fq1_stats = subprocess.check_output("gunzip -c %s | awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'" % (self.fq1_filename, self.fq1_filename), shell=True)
            fq2_stats = subprocess.check_output("gunzip -c %s | awk -v filename=%s '((NR-2)%%4==0){read=$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'" % (self.fq2_filename, self.fq2_filename), shell=True)
        #return subprocess.check_output("parallel \"gunzip -c {} | awk -v filename={} '((NR-2)%%4==0){read=\$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'\" ::: %s" % ' '.join(fq_set))


def format_options(options):
    formatted_options = {}
    for k,v in options.items():
        if k.startswith("-"):
            formatted_options[k] = v
        elif len(k) == 1:
            formatted_options['-' + k] = v
        else:
            formatted_options['--' + k] = v
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
        options = ' '.join([k + " " + v for k,v in format_options(options).items()])
        self.actions += ["bcftools filter %s" % (options)]
        return self

    def stats(self):
        return subprocess.check_output("bcftools stats --samples - %s" % self.filename, stderr=subprocess.STDOUT, shell=True)

    def region(self,chrom,start,end):
        print "bcftools filter -H -r %s:%s-%s %s" % (chrom, start, end, self.file)
        self.actions += ["bcftools view  -r %s:%s-%s %s" % (chrom, start, end, self.file)]
        return self

    def nVariants(self):
        return subprocess.Popen("bcftools view %s | grep -v '^#' | wc -l" % (self.file), shell=True, stdout=subprocess.PIPE).communicate()

    def include(self,depth):
        self.actions += ["bcftools filter --include 'DP<%s'" % (depth)]
        return self

    def out(self, out_filename, type="bcf"):
        # Output types
        actions = ' | '.join(self.actions)
        
        out_command = "bcftools view %s | %s > %s" % (self.filename, actions, out_filename)
        print out_command
        return subprocess.check_output("bcftools view -O b %s | %s > %s" % (self.filename, actions, out_filename))







x = bcf("04_mmp_strains.txt.vcf.gz")
print x.filter({"include":'%QUAL>30', "soft-filter":"MaxQualityFail"})
print x.filter({"include":'DP>3', "soft-filter": "Minimum-Depth"})
print x.actions

x.out("fixed", "bcf")

