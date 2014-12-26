pyPipeline
=========

#### FASTQ

- [ ] Nextera Adapter Trimming
- [ ] Create debug mode, to sample first x number of lines for running analysis.

#### Aligned Sequences: BAM

- [X] Alignment with BWA
- [X] Remove Optical Duplicates (Picard Tools)
- [ ] BAM Statistics Report
- [ ] FASTQC Report
- [ ] Duplicate Reads Report / Summarization

## Variant Calling

- [ ] Joint Calling
	- [ ] Parallelization (by chromosome; then merge)
- [ ] Individual Calling
	- [X] Parallelization (Across Regions)
	- [ ] Merge Individual Samples (Option)
- [ ] Add ability to automatically 'downgrade' VCF for viewing in IGV (temp)

### Callers 

- [ ] Samtools/bcftools
- [ ] Freebayes
- [ ] Platypus

### Filters

- [ ] Soft-Filtering
- [ ] Hard Filtering

### Variant Comparison

- [ ] Depth of Coverage and Mitochondrial Statistics
- [ ] Concordance Analysis and Report
- [ ] Comparing different callers (integrates vcf-toolbox)

### Validation

- [ ] Integrate Swan and Wick variant calls.
- [ ] Integrate original MMP calls.
- [ ] Integrate Radseq calls
- [ ] Develop list of 196 sites to validate; random sample from CB4856.

### Analysis

- [ ] Beagle
- [ ] Deduplication Reports (Picard)
	- [X] Data handling
	- [ ] Report Generation
- [ ] Fastqc Reports
	- [ ] FQ
	- [ ] BAM (individual+merged)

## Other

- [ ] Email Notifications
- [ ] Copy analysis.yaml to analysis_dir
- [ ] Add checks for required options.
- [ ] copy vs. rename single bams.
- [X] Fix gzip vs. bgzip issue.
