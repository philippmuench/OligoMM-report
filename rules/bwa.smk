# mapping QCed reads to joined OMM reference
rule bwa_mem:
	input:
		reads=["processed/trimmed/pe/{sample}.1.fastq.gz", "processed/trimmed/pe/{sample}.2.fastq.gz"]
	output:
		temp("processed/bwa_mapped_omm/{sample}.bam")
	log:
		"logs/bwa_mem/{sample}.log"
	params:
		index=config["REF"],
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="picard",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	threads: 10
	wrapper:
		"0.36.0/bio/bwa/mem"

# mapping QCed PE reads to joined OMM reference + E. Coli
rule bwa_mem_ecoli:
	input:
		reads=["processed/trimmed/pe/{sample}.1.fastq.gz", "processed/trimmed/pe/{sample}.2.fastq.gz"]
	output:
		temp("processed/bwa_mapped_omm_plus_ecoli/{sample}.bam")
	log:
		"logs/bwa_mem_with_ecoli/{sample}.log"
	params:
		index=config["REF_ECOLI"],
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="samtools",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	threads: 
		10
	wrapper:
		"0.36.0/bio/bwa/mem"

rule samtools_index_omm:
	input:
		"processed/bwa_mapped_omm/{sample}.bam"
	output:
		"processed/bwa_mapped_omm/{sample}.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"

rule samtools_index_omm_with_ecoli:
	input:
		"processed/bwa_mapped_omm_plus_ecoli/{sample}.bam"
	output:
		"processed/bwa_mapped_omm_plus_ecoli/{sample}.bam.bai"
	params:
		""
	wrapper:
		"0.36.0/bio/samtools/index"
