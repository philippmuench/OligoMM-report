rule mpilup:
	input:
		# single or list of bam files
		bam="processed/bwa_mapped_omm/{sample}.bam",
		reference_genome=config["REF"]
	output:
		"processed/mpileup/{sample}.mpileup.gz"
	log:
		"logs/samtools/mpileup/{sample}.log"
	wrapper:
		"0.50.4/bio/samtools/mpileup"

rule pileup_to_vcf:
	input:
		gz="processed/mpileup/{sample}.mpileup.gz"
	output:
		snp="processed/varscan/{sample}.snp"
	message:
		"Calling pileup2snp"
	shell:
		"""
		gunzip --stdout {input.gz} | varscan pileup2snp --min-avg-qual 20 --p-value 0.01  > {output.snp}
		"""

rule convertvarscan:
	input:
		snp="processed/varscan/{sample}.snp"
	output:
		vcf="processed/varscan/{sample}.vcf"
	message:
		"Convert Varscan output to vcf"
	log:
		"logs/varscan_{sample}_mpileup2cns.log"
	shell:
		"python ../tools/varscan2vcf.py -s {input.snp} -o {output.vcf}"
