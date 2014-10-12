#! /user/bin/env python
# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, re, glob, hashlib, gzip, mimetypes, tempfile
from itertools import groupby as g
from collections import OrderedDict

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
