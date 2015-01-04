pyPipeline
=========

#### FASTQ

- [ ] Nextera Adapter Trimming
- [X] Create debug mode, to sample first x number of lines for running analysis.

#### Aligned Sequences: BAM

- [X] Alignment with BWA
- [X] Remove Optical Duplicates (Picard Tools)
- [X] Sbatch Job Dependencies
- [ ] Duplicate Reads Report / Summarization

## Variant Calling

- [X] Construct files from fastq list.
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
- [X] Sbatch Job Dependencies
	- [X] Joint Calling
	- [X] Individual Calling
- [ ] Add function for 'downgrading' VCF for viewing in IGV (temporary)

### Callers 

- [X] Samtools/bcftools
- [ ] Freebayes
- [ ] Platypus

### Filters

- [X] Soft-Filtering
- [ ] Hard Filtering

### Validation

- [ ] Integrate Swan and Wick variant calls.
- [ ] Integrate original MMP calls.
- [ ] Integrate Radseq calls
- [ ] Examine TsTv Ratios
- [ ] Develop list of 196 sites to validate; random sample from CB4856.

### Statistics

- [ ] Remove stats files if sample file changed (use hash)
- [X] FASTQ - statistics/
	- [X] Check if stats already taken (based on chksum)
	- [ ] Add date/time
- [ ] BAM - statistics/
	- [ ] Check if stats already taken
	- [ ] Add date/time
- [ ] VCF - statstics/
	- [ ] Check if stats already taken
	- [ ] Add date/time

### Reports

- [ ] Fastq report
- [ ] Deduplication Report
- [ ] Het polarization report
- [ ] BAM report
- [ ] VCF report
- [ ] Variant Comparison
	- [ ] Depth of Coverage and Mitochondrial Statistics
	- [ ] Concordance Analysis and Report
	- [ ] Comparison of different callers

## Other

- [ ] Email Notifications
- [ ] Copy analysis.yaml to analysis_dir
- [ ] Add checks for required options.
- [ ] copy vs. rename single bams.
- [X] Fix gzip vs. bgzip issue.
