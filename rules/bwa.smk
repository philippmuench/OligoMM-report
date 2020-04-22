# mapping QCed reads to joined OMM reference using bwa-mem
rule bwa_mem:
	input:
		reads=["processed/trimmed/pe/{sample}.1.fastq.gz", "processed/trimmed/pe/{sample}.2.fastq.gz"]
	output:
		"processed/{sample}/{sample}_bwa_mapped_omm_sorted.bam"
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

rule index_bwa:
	input:
		"processed/{sample}/{sample}_bwa_mapped_omm_sorted.bam"
	output:
		"processed/{sample}/{sample}_bwa_mapped_omm_sorted.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"
