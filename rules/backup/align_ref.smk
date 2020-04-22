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
	threads: 10
	wrapper:
		"0.36.0/bio/bwa/mem"

# extract unmapped reads
rule extract_unmapped_reads:
	input:
		bam=["mapped/{sample}.bam"]
	output:
		unmapped=["unmapped/{sample}.bam"],
		sorted=["unmapped/{sample}_sorted.bam"]
	shell: """
		samtools view -b -f 4 {input.bam} > {output.unmapped} 
		samtools sort {output.unmapped} -o {output.sorted}
		"""

# bam to fastq to get unmapped PE reads
rule bam_to_fastq:
	input:
		bam=["unmapped/{sample}_sorted.bam"]
	output:
		fq1=["unmapped/fastq/{sample}_1.fastq"],
		fq2=["unmapped/fastq/{sample}_2.fastq"]
	shell: 
		"bedtools bamtofastq -i {input.bam} -fq {output.fq1} -fq2 {output.fq2}"

# assign taxon to non-mapping reads based on kneaddata output
rule kraken2:
	input:
		reads1=["unmapped/fastq/{sample}_1.fastq"],
		reads2=["unmapped/fastq/{sample}_2.fastq"]
	output:
		krak=["kraken/{sample}/{sample}_unmapped.krak"],
		report=["kraken/{sample}/{sample}_unmapped.report"]
	params:
		db=KRAKEN_DB
	threads: 8
	shell: 
		"kraken2 --paired --use-names --db {params.db} --threads 20 --output {output.krak} --report {output.report} {input.reads1} {input.reads2}"

# generate html kraken report
rule kraken2Sankey:
	input:
		report=["kraken/{sample}/{sample}.report"]
	output:
		html=["kraken/{sample}/{sample}.html"]
	shell: 
		"Rscript ../scripts/krakenSankey.R {input.report}"

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
	threads:	# Samtools takes additional threads through its option -@
		25		 # This value - 1 will be sent to -@.
	wrapper:
		"0.36.0/bio/samtools/sort"

rule mark_duplicates:
	input:
		bam="mapped/{sample}.sorted.bam"
	output:
		bam="dedup/{sample}.bam"
	params:
		jar=PICARD_PATH
	log:
		"logs/picard/dedup/{sample}.log"
	shell: 
		"java -Xmx10G -jar {params.jar} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT QUIET=true AS=true METRICS_FILE=metrics VERBOSITY=WARNING TMP_DIR=/scratch/pmuench"

rule samtools_index:
	input:
		"dedup/{sample}.bam"
	output:
		"dedup/{sample}.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"

rule samtools_index_raw:
	input:
		"mapped/{sample}.sorted.bam"
	output:
		"mapped/{sample}.sorted.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.36.0/bio/samtools/index"

# samtools variant calling -------------------------------------------
# use samtool's varaiant calling mpileup
rule mpileup:
	input:
 		bam="mapped_ind/{sample}.{genome}.sorted.bam",
		reference_genome=REF
	output:
		"mpilep_ind/{sample}.{genome}.mpileup.gz"
 	log:
		"logs/samtools/mpileup_ind/{sample}.{genome}.log"
	params:
 		 extra="-d 10000",  # optional
	wrapper:
		"0.39.0/bio/samtools/mpileup"

# convert mpilup output to SNP format
rule mpileup_to_snp:
	input:
		"mpileup_ind/{sample}.{genome}.mpileup.gz"
	output:
		"varscan_ind/{sample}.{genome}.snp"
	log:
		"logs/varscan_ind/{sample}.{genome}.log"
	wrapper:
		"0.39.0/bio/varscan/mpileup2snp"

rule mpileup_to_vcf:
	input:
		mpileup="mpileup_ind/{sample}.{genome}.mpileup.gz"
	output:
		vcf="varscan_ind/{sample}.{genome}.vcf"
	log:
		"logs/varscan_ind/{sample}.{genome}.log"
	shell: "java -Xmx2g -jar ../tools/VarScan.v2.4.4.source.jar mpileup2snp --output-vcf 1 < zcat {input.mpileup} > {output.vcf}"
