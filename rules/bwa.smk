# mapping QCed reads to joined OMM reference
rule bwa_mem:
	input:
		reads=["trimmed/pe/{sample}.1.fastq.gz", "trimmed/pe/{sample}.2.fastq.gz"]
	output:
		temp("mapped/{sample}.bam")
	log:
		"logs/bwa_mem/{sample}.log"
	params:
		index=REF,
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="picard",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	threads: 20
	wrapper:
		"0.36.0/bio/bwa/mem"

# mapping QCed PE reads to joined OMM reference + E. Coli
rule bwa_mem_ecoli:
	input:
		reads=["trimmed/pe/{sample}.1.fastq.gz", "trimmed/pe/{sample}.2.fastq.gz"]
	output:
		temp("mapped_with_e/{sample}.bam")
	log:
		"logs/bwa_mem_with_e/{sample}.log"
	params:
		index=REF_ECOLI,
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="samtools",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	threads: 
		20
	wrapper:
		"0.36.0/bio/bwa/mem"

rule samtools_index:
	input:
		"mapped/{sample}.bam"
	output:
		"mapped/{sample}.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"
