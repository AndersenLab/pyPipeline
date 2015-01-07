import sys, os
from subprocess import Popen, PIPE
from pprint import pprint as pp


class bamfile:
    """
        Basic utility functions
    """

    def __init__(self, filename):
        self.header = self.fetch_header(filename)
        self.RG = self.header["RG"]

    def fetch_header(self, filename):
        comm = ["samtools", "view", "-H", filename]
        out, err = Popen(comm, stdout=PIPE).communicate()
        if err == "":
            raise Exception("Error parsing bam file.")
        else:
            header = out.strip().split("\n")
            header_set = {}
            for tag, val in [x.strip('@').split("\t", 1) for x in header]:
                print tag
                if tag not in header_set:
                    header_set[tag] = []
                val = dict([tuple(x.split(':', 1)) for x in val.split('\t')])
                header_set[tag] += [val]
            header_set[tag] = sorted(header_set[tag])
            return header_set

