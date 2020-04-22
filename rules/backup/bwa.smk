# align PE reads to joined reference
rule bwa_mem:
	input:
		reads=["trimmed/pe/{sample}.1.fastq.gz", "trimmed/pe/{sample}.2.fastq.gz"]
	output:
		"mapped/{sample}.bam"
	log:
		"logs/bwa_mem/{sample}.log"
	params:
		index=REF,
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="none",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	threads: 20
	wrapper:
		"0.36.0/bio/bwa/mem"

rule samtools_flagstat:
	input:
		"mapped/{sample}.bam"
	output:
		"mapped/{sample}.bam.flagstat"
	wrapper:
		"0.36.0/bio/samtools/flagstat"

rule samtools_sort:
	input:
		"mapped/{sample}.bam"
	output:
		"mapped/{sample}.sorted.bam"
	threads:
		10
	wrapper:
		"0.36.0/bio/samtools/sort"

rule mark_duplicates:
	input:
		bam="mapped/{sample}.sorted.bam"
	output:
		bam="mapped_dedup/{sample}.sorted.dedup.bam"
	params:
		jar=PICARD_PATH
	log:
		"logs/picard/dedup/{sample}.log"
	shell: 
		"java -Xmx10G -jar {params.jar} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT QUIET=true AS=true METRICS_FILE=metrics VERBOSITY=WARNING TMP_DIR=/scratch/pmuench"

rule samtools_index:
	input:
		"mapped/{sample}.sorted.bam"
	output:
		"mapped/{sample}.sorted.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"
