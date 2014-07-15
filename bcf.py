# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, uuid, re
import vcf.filters
from collections import OrderedDict


class bcf(file):
    def __init__(self, file):
        # Start by storing basic information about the vcf/bcf
        self.file = file
        self.ops = []
        self.header = subprocess.Popen("bcftools view -h %s" % self.file, shell=True, stdout=subprocess.PIPE).communicate()[0]
        print self.header
        # Samples
        self.samples = filter(len,subprocess.Popen("bcftools query -l %s" % self.file, shell=True, stdout=subprocess.PIPE).communicate()[0].split("\n"))

        # Meta Data
        self.meta = OrderedDict(re.compile(r'''^##(?P<key>[^<#]+?)=(?P<val>[^<#]+$)''', re.M).findall(self.header))

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

    # Parse Header

    def filename(self):
        return self.file

    def region(self,chrom,start,end):
        print "bcftools filter -H -r %s:%s-%s %s" % (chrom, start, end, self.file)
        self.ops += ["bcftools view  -r %s:%s-%s %s" % (chrom, start, end, self.file)]
        return self

    def include(self,depth):
        self.ops += ["bcftools filter --include 'DP<%s'" % (depth)]
        return self

    def out(self):
        print self.ops
        print ' | '.join(self.ops)
        return (len,subprocess.Popen(' | '.join(self.ops + ["bcftools view -H"]), shell=True, stdout=subprocess.PIPE).communicate()[0].split("\n"))







x = bcf("vcf/mmp.vcf.gz")

print x.file
print x.samples
print x.metadata["fileformat"]



#print x.region('chrIII', 3800, 4000).out()