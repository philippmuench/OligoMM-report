# makes a local copy of the read files
rule copy_files:
	input: 
		from1=os.path.join(READS_PATH, "{sample}_R1_001.fastq.gz"),
		from2=os.path.join(READS_PATH, "{sample}_R2_001.fastq.gz")
	output:
		to1="reads/pe/{sample}_R1_001.fastq.gz",
		to2="reads/pe/{sample}_R2_001.fastq.gz"
	shell: """
					cp {input.from1} {output.to1}
					cp {input.from2} {output.to2}
					"""

# QC using fastp (https://github.com/OpenGene/fastp), pipeline will use kneaddata
rule fastp_pe:
		input:
				sample=["reads/pe/{sample}_R1_001.fastq.gz", "reads/pe/{sample}_R2_001.fastq.gz"]
		output:
				trimmed=["trimmed/pe/{sample}.1.fastq.gz", "trimmed/pe/{sample}.2.fastq.gz"],
				html="report/pe/{sample}.html",
				json="report/pe/{sample}.json"
		log:
				"logs/fastp/pe/{sample}.log"
		params:
				extra=""
		threads: 25
		wrapper:
				"0.36.0/bio/fastp"
