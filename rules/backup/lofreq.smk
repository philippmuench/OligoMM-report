# run lofreq on mapping
rule lofreq:
	input:
		bam="mapped/{sample}.sorted.bam",
		bai="mapped/{sample}.sorted.bam.bai",
		ref=REF
	output:
		vcf="lofreq/{sample}.vcf"
	log:
		"logs/lofreq/{sample}.log"
	threads: 40
	shell: 
		"lofreq call-parallel --pp-threads {threads} -f {input.ref} {input.bam} -o {output.vcf} --min-cov 5 --illumina-1.3 --force-overwrite > {log} 2> {log}"

rule bgzip:
	input:
		vcf="lofreq/{sample}.vcf"
	output:
		vcfgz="lofreq/{sample}.vcf.gz"
	log:
		"logs/bgzip/{sample}.log"
	shell:
		"bgzip -c {input.vcf} > {output.vcfgz} 2> {log}"

rule index_vcf:
	input:
		vcfgz="lofreq/{sample}.vcf.gz"
	output:
		vcfsort="lofreq/{sample}.vcf.gz.tbi"
	log:
		"logs/lofreq/index_vcf_{sample}.log"
	shell:
		"tabix -f -p vcf {input.vcfgz}"

rule split_vcf_by_genome:
	input:
		x1="lofreq/{sample}.vcf.gz",
		x2="../databases/oligomm/names/{genome}.txt",
		x3="lofreq/{sample}.vcf.gz.tbi"
	output:
		x3="lofreq_ind/{sample}_{genome}.vcf"
	shell: """
		tabix -h -R {input.x2} {input.x1} > {output.x3}
	"""

#rule get_chr_list:
#	input:
#		vcf="lofreq_ind/{sample}.{genome}.vcf"
#	output:
#		chr="lofreq_ind/{sample}.{genome}.txt"
#	shell:
#		"cat {input.vcf} | mawk '$1 ~ /^#/ {{next}} {{print $1 | "sort -k1,1 -u"}}' > {output.chr}"

#sorts and gz a vcf file and indexes
rule sort_vcf:
	input:
		vcf="lofreq/{sample}.vcf"
	output:
		sorted="lofreq/{sample}.sorted.vcf.gz"
	shell: """
		cat {input.vcf} | mawk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | bgzip -c > {output.sorted}
		tabix -p vcf {output.sorted}
		"""
#rule annotate_vcf:
#	input:
#		vcf="lofreq_ind/{sample}.{genome}.sorted.vcf.gz",
#		bed="../databases/bed/{genome}.bed.gz"
#	output:
#		vcf2="lofreq_annotated/{sample}.{genome}.vcf"
#	shell: 
#		"bcftools annotate -a {input.bed} -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') {input.vcf} > {output.vcf2}"
