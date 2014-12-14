#usr/bin/python

# Align with BWA
bwa = """bwa {bwa_command} -R '{RG_header}' {bwa_options} {reference} {FQ1} {FQ2} | samtools view -bhu - > {OPTIONS.analysis_dir}/bam/{ID}.unsorted.bam
         samtools sort -O bam -T {tmpname} {OPTIONS.analysis_dir}/bam/{ID}.unsorted.bam > {OPTIONS.analysis_dir}/bam/{ID}.sorted.bam"""

mark_dups = """
            java -jar tools/picard.jar MarkDuplicates \
            I={OPTIONS.analysis_dir}/bam/{ID}.sorted.bam \
            O={OPTIONS.analysis_dir}/bam/{ID}.bam \
            M={OPTIONS.analysis_dir}/bam/{ID}.duplicate_report.txt \
            VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
            """

samtools_command = """ """

freebayes_command = """ """
