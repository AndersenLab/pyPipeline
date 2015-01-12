#!/usr/bin/env python
"""pyPipeline

Usage:
  pipe.py trim <config>
  pipe.py align <config> [options]
  pipe.py snps (individual|joint) <config> [options]
<<<<<<< HEAD
  pipe.py transposons <config>
  pipe.py samplefile <filename/dir>
=======
  pipe.py samplefile <config>
>>>>>>> FETCH_HEAD
  pipe.py genome [<name>]
  pipe.py test <config> [options]

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
from utils import *
from utils.constants import *
from utils.configuration import *
from utils.genomes import *
import csv
from pprint import pprint as pp

if __name__ == '__main__':
    opts = docopt(__doc__, version='pyPipeline')
    print opts

    # Setup configuration
    if opts["<config>"] is not None:
        config_file = opts["<config>"]
        cf = config(opts["<config>"])
        # Create new sample file
        if opts["samplefile"] is True:
            if file_exists(cf.sample_file):
                msg("Sample file already exists","error")
            sf = cf.get_sample_file()
            sf.new_sample_file()

    # sample file
    sf = cf.get_sample_file()

    #==================#
    # Genome Retrieval #
    #==================#

    if opts["genome"] is True:
        if opts["<name>"] is not None:
            fetch_genome(opts["<name>"])
        else:
            list_genomes()
        exit(0)

    #=======#
    # Setup #
    #=======#

<<<<<<< HEAD
    analysis_types = ["trim", "align", "merge", "snps", "indels", "test", "transposons"]
=======
>>>>>>> FETCH_HEAD
    analysis_type = [x for x in opts if opts[x] == True and x in analysis_types][0]

<<<<<<< HEAD
    print "{script_dir}/genomes/{OPTIONS.reference}/*fa.gz".format(**locals())
    reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*f*.gz".format(**locals()))[0]
    print reference
    sample_file = open(config.OPTIONS.sample_file, 'rU')

    #======================#
    # Debug (testing) mode #
    #======================#
    if OPTIONS.debug == True:
        debug_fq_dir = "{OPTIONS.fastq_dir}/DEBUG_FQ".format(**locals())
        makedir(debug_fq_dir)
        bam_dir = "{OPTIONS.analysis_dir}/debug".format(**locals())
        eav_file = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}/DEBUG_eav.txt".format(**locals())
=======
    cf.log("#=== Beginning Analysis ===#")
    cf.log("Running " + opts["<config>"])
>>>>>>> FETCH_HEAD

    #======#
    # Trim #
    #======#

    if analysis_type == "trim":
        # Trim Nextera Adapters
        pass

    #===========#
    # Alignment #
    #===========#

    elif analysis_type == "align":
        cf.log("Performing Alignment")
        for bam in sf.check_bams():
            # Check merged bam
            dependency_list = []
            if bam["bam_merged_exists_and_RG_correct"] is False:
                # Remove merged bam if it exists (RG is wrong)
                if file_exists(bam["bam_merged_filename"]):
                    cf.log("merged bam %s exists or incorrect read group." %
                    bam["bam_merged_filename"])
                    cf.command("rm %s" % bam["bam_merged_filename"])
                # Check individual bams
                for ind_bam, ind_bam_exists, fq, RG, ID in zip(bam["bam_ind_filename"],
                                                                   bam["bam_ind_exists_and_RG_correct"],
                                                                   bam["fq"],
                                                                   bam["raw_RG"],
                                                                   bam["ID"]):
                    if ind_bam_exists is False:
                        fq = {"fq1": fq[0], "fq2": fq[1], "ID": ID, "RG": RG, "SM": bam["SM"]}
                        align = "{run} {script_dir}/align.py {config_file} \"{fq}\"".format(**locals())
                        jobid = cf.submit_job(align,
                                              analysis_type=analysis_type,
                                              log_name=ID)
                        dependency_list.append(jobid)
                # Merge Bams
                bams_to_merge = bam["bam_ind_filename"]
                merge_bams = "{run} {script_dir}/merge_bams.py {config_file} \"{bams_to_merge}\"".format(**locals())
                print pp(bam)
                jobid = cf.submit_job(merge_bams,
                                      log_name=bam["SM"],
                                      analysis_type=analysis_type,
                                      dependencies=dependency_list,
                                      dependency_type="afterok")

    #=============#
    # SNP Calling #
    #=============#
    elif analysis_type == "snps" and opts["individual"] == True:
        #
        # Individual - Needs to be run twice to generate union output
        #

        # Get list of bams
        log.info("Performing Variant Calling")

        # Construct Sample Set
        bam_set = []
        for row in csv.DictReader(sample_file, delimiter='\t', quoting=csv.QUOTE_NONE):
            SM = row["SM"]
            bam_set.append(SM)
            bam_file = "{bam_dir}/{SM}.bam".format(**locals())
            if (not file_exists(bam_file) or file_exists(bam_file + ".csi")) and OPTIONS.debug == False:
                raise Exception("Bam File or index does not exist: %s" % bam_file)
        
        bam_set = set(bam_set)
        # Has vcf been called for given snp caller?
        dependency_list = []
        for SM in bam_set:
            for caller in snp_callers:
                union_vcf_file = "{vcf_dir}/{SM}.{caller}.union.vcf.gz".format(**locals())
                union_vcf_index_file = union_vcf_file + ".csi"
                variant_files = [union_vcf_file, union_vcf_index_file]
                union_variant_file = "{OPTIONS.analysis_dir}/{caller}.{OPTIONS.union_variants}.txt".format(**locals())
                complete_individual = "{vcf_dir}/{SM}.bcftools.individual.vcf.gz".format(**locals())
                individual_vcfs_check = (not file_exists(complete_individual) and COMMANDS.snps.snp_options.remove_temp == False)
                if not all(map(file_exists, variant_files )) or not file_exists(union_variant_file) or individual_vcfs_check:
                    call_snps = """{run} {output_dirs} {script_dir}/call_snps_individual.py {config_file} \"{SM}.bam\"""".format(**locals())
                    jobid = submit_job(call_snps)
                    dependency_list.append(jobid)
        # Merge individual bams
        if COMMANDS.snps.snp_options.merge_individual_vcfs == True:
            merge_snps = """{run} {script_dir}/merge_vcfs_individual.py {config_file}""".format(**locals())
            submit_job(merge_snps, dependency_list)

    elif analysis_type == "snps" and opts["joint"] == True:
        #
        # Joint
        #
        dependency_list = []
        chunks = [x for x in chunk_genome(OPTIONS.chrom_chunk_kb,reference)]
        for caller in snp_callers:
            merged_vcf_name = "{vcf_dir}/joint.{caller}.vcf.gz".format(**locals())
            merged_file_exists = file_exists(merged_vcf_name) and file_exists(merged_vcf_name + ".csi")
            for chunk in chunks:
                # Check that chunk does not exist.
                chunk_sanitized = chunk.replace(":","_")
                vcf_file = "{vcf_dir}/TMP.joint.{chunk_sanitized}.{caller}.vcf.gz".format(**locals())
                if (not file_exists(vcf_file) or not file_exists(vcf_file + ".csi")) and not merged_file_exists:
                    call_snps = """{run} {script_dir}/call_snps_joint.py {config_file} \"{chunk}\"""".format(**locals())
                    jobid = submit_job(call_snps)
                    dependency_list.append(jobid)
        if not merged_file_exists:
            merge_snps = """{run} {script_dir}/concat_vcfs_joint.py {config_file}""".format(**locals())
            submit_job(merge_snps, dependency_list)
        else:
            print("Merged File Already Exists")
    elif analysis_type == "test":
        print "GREAT"
<<<<<<< HEAD
        r = "{run} {script_dir}/call_snps_individual.py {config_file} '[1,2,3]'".format(**locals())
        os.system(r)
    elif analysis_type == "transposons":
        bam_set = []
        for fq in csv.DictReader(sample_file, delimiter='\t', quoting=csv.QUOTE_NONE):
            if fq["RUN"] != "NO":
                SM = fq["SM"]
                bam_set.append(SM)
                bam_file = "{bam_dir}/{SM}.bam".format(**locals())
        bam_set = set(bam_set)
        # Has vcf been called for given snp caller?
        dependency_list = []
        for SM in bam_set:
            call_transposons = """{run} {script_dir}/call_transposons.py {config_file} \"{SM}.bam\"""".format(**locals())
            print call_transposons
            jobid = submit_job(call_transposons)
            #dependency_list.append(jobid)



=======
>>>>>>> FETCH_HEAD


