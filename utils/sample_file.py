




class sample_file:
	"""
		Class for handling actions associated with the sample file.
	"""
	def __init__(self, filename):
		self.sample_file = open(filename, 'r')

		# Run qc on sample file
		for fq in csv.DictReader(sample_file, delimiter='\t', quoting=csv.QUOTE_NONE)
			if fq["RUN"] != "NO":
                fq1, fq2 = fq["FQ1"], fq["FQ2"]
                fq["fq1"] = "{OPTIONS.fastq_dir}/{fq1}".format(**locals())
                fq["fq2"] = "{OPTIONS.fastq_dir}/{fq2}".format(**locals())
                # Construct Individual BAM Dict
                ID = fq["ID"]
                SM = fq["SM"]

                # Sanity Checks
                if ID in ID_set:
                    raise Exception("IDs are not unique: %s" % ID)
                else:
                    ID_set.append(ID)
                if SM not in sample_set:
                    sample_set[SM] = []
                if ID == None or ID == "":
                    raise Exception("No ID defined for %s" % fq)
                if SM == None or SM == "":
                    raise Exception("No sample defined for %s" % ID)
                if fq1 == fq2:
                    raise Exception("Both Fastq's share same name: %s %s" % (fq1, fq2))
                if ID == SM:
                    raise Exception("ID cannot be equal to SM")

                RG = construct_RG_header(ID, fq).replace("\\t","\t")
                sample_info = {"ID" : ID, "RG": RG, "fq": fq}
                sample_set[SM].append(sample_info)

        def process

s = sample_file("/Users/dancook/Documents/tmp/data/FASTQ.txt")

print s.sample_file.read()