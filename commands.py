#usr/bin/python

# Align with BWA
bwa = """bwa {bwa_command} -R '{RG_header}' {bwa_options} {reference} {FQ1} {FQ2} | samtools view -bhu - > {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.unsorted.bam
         samtools sort -O bam -T {tmpname} {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.unsorted.bam > {OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam"""

mark_dups = """
            java -jar {script_dir}/tools/picard.jar MarkDuplicates \
            I={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.sorted.bam \
            O={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.bam \
            M={OPTIONS.analysis_dir}/{OPTIONS.bam_dir}/{ID}.duplicate_report.txt \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
            """

merge_bams = """samtools merge {merge_options} {merged_bam_name} {constitutive_bams}"""

samtools_command = """ """

freebayes_command = """ """
