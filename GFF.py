#! /user/bin/env python
# wrapper for tabix around bcbio gff parser; allows one to query a gff file.
import os, subprocess
from BCBio.GFF.GFFParser import _gff_line_map

class dummy_params(object):
    def __init__(self):
        self.limit_info = {}
        self.jsonify = False

class gff(object):
    def __init__(self, filename):
        # Fetch the gff header
        self.header = subprocess.check_output(["tabix","-h",filename,"chr:0-0"]).strip().split("\n")
        self.filename = filename   

    def query(self, chrom, start, end):
        """
        Call tabix and generate an array of strings for each line it returns.
        * * * Thanks to https://github.com/slowkow/pytabix * * *
        """
        query = '{}:{}-{}'.format(chrom, start, end)
        process = subprocess.Popen(['tabix', '-f', self.filename, query], stdout=subprocess.PIPE)
        for line in process.stdout:
            yield _gff_line_map(line.strip(), dummy_params())