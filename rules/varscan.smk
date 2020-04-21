rule mpilup:
	input:
		# single or list of bam files
		bam="mapped/{sample}.bam",
		reference_genome=REF
	output:
		"mpileup/{sample}.mpileup.gz"
	log:
		"logs/samtools/mpileup/{sample}.log"
	wrapper:
		"0.50.4/bio/samtools/mpileup"

rule pileup_to_vcf:
	input:
		gz="mpileup/{sample}.mpileup.gz"
	output:
		snp="vcf_varscan/{sample}.snp"
	message:
		"Calling pileup2snp"
	shell:
		"""
		gunzip --stdout {input.gz} | varscan pileup2snp --min-avg-qual 20 --p-value 0.01  > {output.snp}
		"""

#rule mpileup_to_vcf:
#	input:
#		"mpileup/{sample}.mpileup.gz"
#	output:
#		"vcf_varscan/{sample}.snp"
#	message:
#		"Calling SNP with Varscan2"
#	threads:  # Varscan does not take any threading information
#		1     # However, mpileup might have to be unzipped.
#	log:
#		"logs/varscan_{sample}.log"
#	wrapper:
#		"0.50.4/bio/varscan/mpileup2snp"

rule convertvarscan:
	input:
		snp="vcf_varscan/{sample}.snp"
	output:
		vcf="vcf_varscan/{sample}.vcf"
	message:
		"Convert Varscan output to vcf"
	log:
		"logs/varscan_{sample}_mpileup2cns.log"
	shell:
		"python ../tools/varscan2vcf.py -s {input.snp} -o {output.vcf}"
