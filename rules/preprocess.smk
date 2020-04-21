# makes a local copy of the reads
rule copy_files:
	input: 
		from1=os.path.join(READS_PATH, "{sample}_R1_001.fastq.gz"),
		from2=os.path.join(READS_PATH, "{sample}_R2_001.fastq.gz")
	output:
		to1=temp("processed/reads/pe/{sample}_R1_001.fastq.gz"),
		to2=temp("processed/reads/pe/{sample}_R2_001.fastq.gz")
	conda:
		"../envs/conda.yaml"
	shell: """
					cp {input.from1} {output.to1}
					cp {input.from2} {output.to2}
					"""

# QC using fastp (https://github.com/OpenGene/fastp)
rule fastp:
	input:
			sample=["processed/reads/pe/{sample}_R1_001.fastq.gz", "processed/reads/pe/{sample}_R2_001.fastq.gz"]
	output:
		trimmed=[temp("processed/trimmed/pe/{sample}.1.fastq.gz"), temp("processed/trimmed/pe/{sample}.2.fastq.gz")],
		html="docs/reports/fastp/{sample}.html",
		json="processed/fastp/{sample}.json"
	conda:
		"../envs/conda.yaml"
	log:
		"logs/fastp/pe/{sample}.log"
	params:
			extra=""
	threads: 25
	wrapper:
		"0.36.0/bio/fastp"
