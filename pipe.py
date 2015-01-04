#!/usr/bin/env python
"""pyPipeline

Usage:
  pipe.py trim <config>
  pipe.py align <config> [options]
  pipe.py snps (individual|joint) <config> [options]
  pipe.py samplefile <filename/dir>
  pipe.py genome [<name>]
  pipe.py test <config> [options]

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
import glob
from utils import *
from utils.genomes import *
import csv
import operator
from pprint import pprint as pp

def check_fqs(fq):
    if not file_exists(fq["fq1"]) or not file_exists(fq["fq2"]):
        raise Exception("File Missing; Check: {fq1}, {fq2}".format(fq1=fq["fq1"], fq2=fq["fq2"]))

def submit_job(command, dependencies = None, dep_type = "afterok"):
    log.info(command)
    if LOCAL == False:
        if dependencies is not None:
            if len(dependencies) > 0:
                dependencies = ':'.join(dependencies)
                depends_on = "--dependency={dep_type}:".format(**locals())
                depends_on += dependencies
            else:
                depends_on = ""
            command = command.split(" ")
            command.insert(1, depends_on)
            command = ' '.join(command)
        else:
            depends_on = ""
        print command
        jobid, err = Popen(command, stdout=PIPE, stderr=PIPE, shell=True).communicate()
        jobid = jobid.strip().split(" ")[-1]
        if jobid.isdigit() == False:
            raise Exception("Error submitting %s" % jobid)
            exit()
        else:
            return jobid
    else:
        os.system(command)
        return None

if __name__ == '__main__':
    opts = docopt(__doc__, version='pyPipeline')

    #========#
    # Basics #
    #========#

    config_file = opts["<config>"]
    analysis_dir = os.getcwd()

    #
    # Add Checks here for required options
    #

    #==================#
    # Genome Retrieval #
    #==================#

    if opts["genome"] == True:
        if opts["<name>"] is not None:
            fetch_genome(opts["<name>"])
        else:
            list_genomes()
        exit()

    #=================#
    # New Sample File #
    #=================#

    """
        Create a sample file where Library, Sample, and Platform info can be added.
        Optionally update an analysis file
    """

    if opts["samplefile"] == True: 
        header = "FQ1\tFQ2\tID\tLB\tSM\tPL\tRUN\n"
        new_sample_file = open(opts["<filename/dir>"] + ".txt",'w')
        new_sample_file.write(header)
        if is_dir(opts["<filename/dir>"]):
            # Construct a sample file using the directory info.
            sample_set = sorted(glob.glob(opts["<filename/dir>"] + "/*.fq.gz"))
            fastq_pairs = zip(sorted([os.path.split(x)[1] for x in sample_set if x.find("1.fq.gz") != -1]), \
                sorted([os.path.split(x)[1] for x in sample_set if x.find("2.fq.gz") != -1]))
            for pair in fastq_pairs:
                ID = get_fq_ID(pair)
                new_sample_file.write("\t".join(pair) + "\t" + ID + "\n")
        exit()

    #=======#
    # Setup #
    #=======#

    analysis_types = ["trim", "align", "merge", "snps", "indels", "test"]
    analysis_type = [x for x in opts if opts[x] == True and x in analysis_types][0]
    # Load Configuration
    config, log, c_log = load_config_and_log(config_file, analysis_type)
    OPTIONS = config.OPTIONS
    COMMANDS = config.COMMANDS

    # Running locally or on a cluster
    log_dir = "{OPTIONS.analysis_dir}/log".format(**locals())
    stat_dir = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}".format(**locals())
    
    makedir(OPTIONS.analysis_dir)
    makedir(log_dir)
    makedir(stat_dir)

    if LOCAL == True:
        run = "python"
        log_files = ""
    else:
        run = "sbatch "

    bam_dir = "{OPTIONS.analysis_dir}/{OPTIONS.bam_dir}".format(**locals())
    vcf_dir = "{OPTIONS.analysis_dir}/{OPTIONS.vcf_dir}".format(**locals())
    
    snp_callers = COMMANDS.snps
    snp_callers = [x for x in snp_callers if x in available_snp_callers]

    log.info("#=== Beginning Analysis ===#")
    log.info("Running " + opts["<config>"])

    reference = glob.glob("{script_dir}/genomes/{OPTIONS.reference}/*fa.gz".format(**locals()))[0]
    sample_file = open(config.OPTIONS.sample_file, 'rU')

    #======================#
    # Debug (testing) mode #
    #======================#
    if OPTIONS.debug == True:
        debug_fq_dir = "{OPTIONS.fastq_dir}/DEBUG_FQ".format(**locals())
        makedir(debug_fq_dir)
        bam_dir = "{OPTIONS.analysis_dir}/debug".format(**locals())
        eav_file = "{OPTIONS.analysis_dir}/{OPTIONS.stat_dir}/DEBUG_eav.txt".format(**locals())

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
        log.info("Performing Alignment")
        sample_set = {} # Generate a list of samples.
        bam_dir_white_list = [] # List of bams to keep following alignment; removes extras
        # Construct Sample Set
        ID_set = [] # Used to check uniqueness of IDs
        for fq in csv.DictReader(sample_file, delimiter='\t', quoting=csv.QUOTE_NONE):
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

                RG = construct_RG_header(ID, fq).replace("\\t","\t")
                sample_info = {"ID" : ID, "RG": RG, "fq": fq}
                sample_set[SM].append(sample_info)

                # Check that fq's exist before proceeding.
                if OPTIONS.debug == False:
                    check_fqs(fq)

        if len(ID_set) > len(set(ID_set)):
            raise Exception("ID's are not Unique")

        dependency_list = {} # Used to keep jobs working in the proper order.
        for SM in sample_set.keys():
            dependency_list[SM] = []
            # Check the header of the merged bam to see if 
            # current file already exists within
            completed_merged_bam = "{bam_dir}/{SM}.bam".format(**locals())
            # Check to see if merged bam contains constitutive bams
            if file_exists(completed_merged_bam) and file_exists(completed_merged_bam + ".bai"):
                RG = get_bam_RG(completed_merged_bam)

                RG_ind = [x["RG"] for x in sample_set[SM]]
                if split_read_group(RG_ind) != split_read_group(RG):
                    # Delete merged Bam, and re-align all individual.
                    log.info("RG do not match; deleting.")
                    remove_file(completed_merged_bam)
                else:
                    log.info("{SM}.bam contains all specified individual bams.".format(**locals()))

            # Align fastq sets
            if not file_exists(completed_merged_bam) or not file_exists(completed_merged_bam + ".bai"):
                for seq_run in sample_set[SM]:
                    ID = seq_run["ID"]
                    fq = seq_run["fq"]
                    single_bam = "{bam_dir}/{ID}.bam".format(**locals())
                    if OPTIONS.alignment_options.remove_temp == False:
                        bam_dir_white_list.append(single_bam)
                        bam_dir_white_list.append(single_bam + ".bai")
                    # Check single bam RG
                    re_align = False
                    if file_exists(single_bam):
                        current_RG = sorted(get_bam_RG(single_bam)[0].split("\t"))
                        seq_run["RG"] = sorted(seq_run["RG"].split("\t"))
                        single_RG_incorrect = (seq_run["RG"] != current_RG)
                        if (seq_run["RG"] != current_RG):
                            log.info("Readgroup for {single_bam} does not match file; deleting".format(**locals()))
                            remove_file(single_bam)
                            remove_file(single_bam + ".bai")
                            re_align = True

                    if not file_exists(single_bam) or re_align:

                        if LOCAL == False:
                            run += " --output={log_dir}/align.{ID}.%j.txt --error={log_dir}/align.{ID}.%j.err ".format(**locals())

                        align = "{run} {script_dir}/align.py {config_file} \"{fq}\"".format(**locals())
                        jobid = submit_job(align)
                        dependency_list[SM].append(jobid)
                    else:
                        log.info("%-50s already aligned individually, skipping" % single_bam)

        #
        # Merging
        #
        for SM in sample_set.keys():  
            completed_merged_bam = "{bam_dir}/{SM}.bam".format(**locals())
            bam_dir_white_list.append(completed_merged_bam)
            bam_dir_white_list.append(completed_merged_bam + ".bai")
            # Merge Bams for same samples.
            if not file_exists(completed_merged_bam) or not file_exists(completed_merged_bam + ".bai"):
                bam_set = [x["ID"] + ".bam" for x in sample_set[SM]]
                bams_to_merge = (SM, bam_set)

                merge_bams = "{run} {script_dir}/merge_bams.py {config_file} \"{bams_to_merge}\"".format(**locals())
                jobid = submit_job(merge_bams, dependency_list[SM], "afterok")
                print jobid
                if LOCAL == False:
                    print("Submitted merge:{SM}; depends on: {ls}".format(SM=SM, ls=','.join(dependency_list[SM])))
            else:
                log.info("%-50s already exists with all individual bams, skipping" % completed_merged_bam)

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
        for fq in csv.DictReader(sample_file, delimiter='\t', quoting=csv.QUOTE_NONE):
            if fq["RUN"] != "NO":
                SM = fq["SM"]
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
                    call_snps = """{run} {script_dir}/call_snps_individual.py {config_file} \"{SM}.bam\"""".format(**locals())
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

    else:
        exit()
    if analysis_type == "test":
        print "GREAT"
        r = "{run} {script_dir}/call_snps_individual.py {config_file} '[1,2,3]'".format(**locals())
        os.system(r)
        





    
