import sys,os
from subprocess import Popen, PIPE
from pprint import pprint as pp

class bamfile:
	"""
		Basic Bam Functions
	"""

	def __init__(self, filename):
		self.header = self.fetch_header(filename)
		self.RG = 1

	def fetch_header(self, filename):
		out, err = Popen(["samtools","view","-H",filename], stdout=PIPE).communicate()
		if err == "":
			raise Exception("Error parsing bam file.")
		else:
			header = out.strip().split("\n")
			header_set = {}
			for tag, val in [x.strip('@').split("\t",1) for x in header]:
				print tag
				if tag not in header_set:
					header_set[tag] = []
				val = dict([tuple(x.split(':',1)) for x in val.split('\t')])
				header_set[tag] += [val]
			header_set[tag] = sorted(header_set[tag])
			return header_set



x = bamfile("/Users/dancook/Documents/tmp/data/testing2/bam/test1.bam")

print pp(x.header)