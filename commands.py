#usr/bin/python

# Align with BWA

mark_dups = """
            java -jar {script_dir}/tools/picard.jar MarkDuplicates \
            I={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam \
            O={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam \
            M={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.duplicate_report.txt \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
            """

merge_bams = """samtools merge {merge_options} {bam_dir}/{merged_bam_name} {bam_dir}/{SM_Bams}"""

samtools_command = """ """

freebayes_command = """ """
