# mapping QCed reads to joined OMM reference using bowtie2
rule bowtie2:
	input:
		sample=["processed/trimmed/pe/{sample}.1.fastq.gz", "processed/trimmed/pe/{sample}.2.fastq.gz"]
	output:
		temp("processed/{sample}/{sample}_bowtie2_mapped_omm.bam")
	log:
		"logs/bowtie2/{sample}.log"
	params:
		index="/net/sgi/oligomm_ab/OligoMM-report/databases/omm/joined_reference",  # prefix of reference genome index (built with bowtie2-build)
 		extra=""  # optional parameters
	threads: 10  # Use at least two threads
	wrapper:
		"0.51.3/bio/bowtie2/align"

rule samtools_sort_bowtie2:
	input:
		"processed/{sample}/{sample}_bowtie2_mapped_omm.bam"
	output:
		"processed/{sample}/{sample}_bowtie2_mapped_omm_sorted.bam"
	params:
		"-m 10G"
	threads:  # Samtools takes additional threads through its option -@
		8
	wrapper:
		"0.51.3/bio/samtools/sort"

# index
rule samtools_index_bowtie:
	input:
		"processed/{sample}/{sample}_bowtie2_mapped_omm_sorted.bam"
	output:
		"processed/{sample}/{sample}_bowtie2_mapped_omm_sorted.bam.bai"
	params:
		""
	wrapper:
		"0.36.0/bio/samtools/index"

