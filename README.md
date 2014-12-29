pyPipeline
=========

#### FASTQ

- [ ] Nextera Adapter Trimming
- [ ] Create debug mode, to sample first x number of lines for running analysis.

#### Aligned Sequences: BAM

- [X] Alignment with BWA
- [X] Remove Optical Duplicates (Picard Tools)
- [X] Sbatch Job Dependencies
- [ ] BAM Statistics Report
- [ ] FASTQC Report
- [ ] Duplicate Reads Report / Summarization

## Variant Calling

- [ ] Construct files from fastq list.
- [X] Joint Calling
	- [X] Parallelization (by chromosome; then merge)
	- [ ] Filters
		- [X] Soft
		- [ ] Hard
	- [X] Heterozygous Polarization Filter
- [X] Individual Calling
	- [X] Parallelization (Across Regions)
	- [X] Merge Individual Samples (Option)
	- [X] Standard Filters
		- [ ] Hard Filters
	- [X] Remove individual(temp files; option)
	- [X] Heterozygous Polarization Filter
	- [X] Generate joint sites
	- [X] Recall with joint variant set
- [ ] Sbatch Job Dependencies
- [ ] Add function for 'downgrading' VCF for viewing in IGV (temporary)

### Callers 

- [X] Samtools/bcftools
- [ ] Freebayes
- [ ] Platypus

### Filters

- [X] Soft-Filtering
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
	- [ ] Sort and make unique during analysis phase.
	- [ ] Report Generation
- [ ] Heterozygous Polarization Reports
- [ ] Fastqc Reports
	- [ ] FQ
	- [ ] BAM (individual+merged)
- [ ] Comparative Reports

## Other

- [ ] Email Notifications
- [ ] Copy analysis.yaml to analysis_dir
- [ ] Add checks for required options.
- [ ] copy vs. rename single bams.
- [X] Fix gzip vs. bgzip issue.
