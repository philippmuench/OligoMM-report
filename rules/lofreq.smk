# run lofreq on mapping
rule lofreq_bwa:
	input:
		bam="processed/bwa_mapped_omm/{sample}.bam",
		bai="processed/bwa_mapped_omm/{sample}.bam.bai",
		ref=config["REF"]
	output:
		vcf=temp("processed/lofreq_bwa/{sample}.vcf")
	log:
		"logs/lofreq_bwa/{sample}.log"
	threads: 20
	shell: 
		"lofreq call-parallel --pp-threads {threads} -f {input.ref} {input.bam} -o {output.vcf} --min-cov 5 --illumina-1.3 --force-overwrite > {log} 2> {log}"

rule index_vcf_bwa:
	input:
		vcfgz="processed/lofreq_bwa/{sample}.vcf.gz"
	output:
		vcfsort="processed/lofreq_bwa/{sample}.vcf.gz.tbi"
	log:
		"logs/lofreq_bwa/index_vcf_{sample}.log"
	shell:
		"tabix -f -p vcf {input.vcfgz}"

#irule split_vcf_by_genome:
#	input:
#		x1="processed/lofreq/{sample}.vcf.gz",
#		x2="../databases/oligomm/names/{genome}.txt",
#		x3="lofreq/{sample}.vcf.gz.tbi"
#	output:
#		x3="lofreq_ind/{sample}_{genome}.vcf"
#	shell: """
#		tabix -h -R {input.x2} {input.x1} > {output.x3}
#	"""

#sorts and gz a vcf file and indexes
rule sort_vcf_bwa:
	input:
		vcf="processed/lofreq_bwa/{sample}.vcf"
	output:
		sorted="processed/lofreq_bwa_sorted/{sample}.vcf.gz"
	shell: """
		cat {input.vcf} | mawk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | bgzip -c > {output.sorted}
		tabix -p vcf {output.sorted}
		"""
