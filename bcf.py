#! /user/bin/env python
# Small wrapper for bcftools 2.0 - for querying vcf data quickly and easily.
import os, subprocess, re, glob, hashlib, gzip, mimetypes, tempfile
from itertools import groupby as g
from collections import OrderedDict
from utils import *
from pprint import pprint as pp

# Variable types defined in vcf header
_vcf_variable_types = {"Integer" : int, "String": str, "Float" : float, "Flag": bool}

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


class vcf(file):
    def __init__(self, filename, region = ""):
        # Start by storing basic information about the vcf/bcf and checking that an index exists.
        self.filename = filename
        self.filetype = os.path.splitext(filename)[1]
        self.region = region
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
        r = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.info_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}
        
        # Filter
        r = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}
        
        r = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.format_set = {x["id"]:x for x in [m.groupdict() for m in r.finditer(self.header)]}

    def parse_info(self, info_line):
        info = [x.split("=") for x in info_line.split(";")]
        info_line_dict = {}
        # Set variable types if available
        for i in info:
            if i[0] in self.info_set:
                # Check if the field consists of multiple values
                if int(self.info_set[i[0]]["number"]) > 1:
                    info_line_dict[i[0]] = map(_vcf_variable_types[self.info_set[i[0]]["type"]], i[1].split(","))
                # Check if field has a value (non-bool)
                elif len(i) > 1:
                    info_line_dict[i[0]] = map(_vcf_variable_types[self.info_set[i[0]]["type"]], [i[1]])[0]
                else:
                    info_line_dict[i[0]] = i[0]
        return info_line_dict

    def parse_format(self, format_set):
        geno_dict = {}
        format = format_set[0].split(":")
        for k,rec in enumerate(format_set[1:]):
            sample_geno = dict(zip(format,rec.split(":")))
            
            # Set Variable Types
            for k,v in sample_geno.items():
                print k, v
                # Check if the field consists of multiple values
                if int(self.format_set[k]["number"]) > 1:
                    sample_geno[k] = map(_vcf_variable_types[self.format_set[k]["type"]], rec.split(","))
                # Check if field has a value (non-bool)
                elif int(self.format_set[k]["number"]) == 1:
                    sample_geno[k] = map(_vcf_variable_types[self.format_set[k]["type"]], [rec])[0]
                else:
                    sample_geno[k] = k



            sample_geno["sample"] = self.samples[k]

    def variant_dict(self, lines=10):
        """
            Parses a set of variants and inserts results into a dictionary.
        """
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
            variant["QUAL"] = float(line[5])
            variant["FILTER"] = line[6].split(";")
            # Parse INFO and set types
            variant["INFO"] = self.parse_info(line[7])
            #variant["GENO"] = self.parse_format(line[8:])
            variant_set.append(variant)
        print pp(variant_set)
        return variant_set

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
        self.region = region
        return self

    def snpeff(self, annotation_db):
        # Apply snpeff annotations
        self.actions += ["bcftools view | snpeff eff %s" % annotation_db]
        self.header_add_lines += ["##INFO=<ID=EFF,Number=1,Type=String,Description=\"SNPEFF Annotation\">"]
        return self

    def rename_samples(self, sample_names):
        # Rename Samples by specifying a dictionary.
        self.samples = [sample_names[x] if x in sample_names else x for x in self.samples]
        return self

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
        self.header[-1] = '\t'.join(self.header[-1].split("\t")[:9] + self.samples) 
        self.header = '\n'.join(self.header).strip()
        tmp_header.write(self.header)
        tmp_header.flush() # Flush out new header

        #====================================#
        # Add filters, reheader, annotations #
        #====================================#

        self.actions += ["bcftools reheader -h %s" % tmp_header.name]


        if version == 4.1:
            self.actions += ["bcftools view | grep -v '##INFO' | bcftools view -O v"]


        # Determine Output file type
        out_ext = os.path.splitext(out_filename)[1]
        out_opts = { ".bcf" : "b", ".gz" : "z", ".vcf" : "v"}
        if out_ext in out_opts.keys():
            self.actions += ["bcftools view -O %s" % out_opts[out_ext]]
        else:
            raise Exception("Unknown file extension (%s); Must Specify vcf, vcf.gz, or bcf" % out_filename)


        # Output types
        actions = ' | '.join(self.actions)

        print "bcftools view -O u %s %s | %s > %s" % (self.filename, self.region, actions, out_filename)
        subprocess.check_output("bcftools view -O u %s | %s > %s" % (self.filename, actions, out_filename), shell=True)
        if self.filetype in [".bcf", ".gz"]:
            subprocess.check_output("bcftools index %s" % out_filename, shell=True)


#x = bcf("NIC276.nofilter.group.bcf")
#x.filter({"include":'%QUAL>30', "soft-filter":"MaxQualityFail"})
#x.filter({"include":'DP>3', "s": "MinimumDepth"})
#x.out("NIC276.vcf.gz", version = 4.2)

