# run lofreq on mapping to combined reference
rule lofreq:
	input:
		bam="mapped/{sample}.sorted.bam",
		bai="mapped/{sample}.sorted.bam.bai",
		ref=REF
	output:
		vcf="lofreq_call/{sample}.vcf"
	log:
		"logs/lofreq_call/{sample}.log"
	threads: 20
	shell: 
		"lofreq call-parallel --pp-threads {threads} -f {input.ref} {input.bam} -o {output.vcf} > {log} 2> {log}"

# use samtool's varaiant calling mpileup
rule mpileup:
	input:
 		bam="mapped/{sample}.sorted.bam",
		reference_genome=REF
	output:
		"mpilep/{sample}.mpileup.gz"
 	log:
		"logs/samtools/mpileup/{sample}.log"
	params:
 		 extra="-d 10000",  # optional
	wrapper:
		"0.39.0/bio/samtools/mpileup"

# convert mpilup output to SNP format
rule mpileup_to_snp:
	input:
		"mpileup/{sample}.mpileup.gz"
	output:
		"varscan/{sample}.snp"
	log:
		"logs/varscan_{sample}.log"
	wrapper:
		"0.39.0/bio/varscan/mpileup2snp"

rule mpileup_to_vcf:
	input:
		mpileup="mpileup/{sample}.mpileup.gz"
	output:
		vcf="varscan/{sample}.vcf"
	log:
		"logs/varscan_vcf_{sample}.log"
	shell: "java -Xmx2g -jar ../tools/VarScan.v2.4.4.source.jar mpileup2snp --output-vcf 1 < zcat {input.mpileup} > {output.vcf}"

rule bgzip:
	input:
		vcf="lofreq_call/{sample}.vcf"
	output:
		vcfgz="lofreq_call/{sample}.vcf.gz"
	log:
		"logs/bgzip/{sample}.log"
	shell:
		"bgzip -c {input.vcf} > {output.vcfgz} 2> {log}"

rule index_vcf:
	input:
		vcfgz="lofreq_call/{sample}.vcf.gz"
	output:
		vcfsort="lofreq_call/{sample}.vcf.gz.tbi"
	log:
		"logs/sort_vcf/{sample}.log"
	shell:
		"tabix -f -p vcf {input.vcfgz}"

rule get_chr_list:
	input:
		vcf="lofreq_call/{sample}.vcf"
	output:
		chr="lofreq_call/{sample}.txt"
	shell:
		"cat {input.vcf} | mawk '$1 ~ /^#/ {{next}} {{print $1 | "sort -k1,1 -u"}}' > {output.chr}"

#sorts and gz a vcf file and indexes

rule sort_vcf:
	input:
		vcf="lofreq_call/{sample}.vcf"
	output:
		sorted="lofreq_call/{sample}_sorted.vcf.gz"
	shell: """
		cat {input.vcf} | mawk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | bgzip -c > {output.sorted}
		tabix -p vcf {output.sorted}
		"""

#extracts a genomes from a sorted gz vcfDD

rule split_vcf:
	input:
		vcf="lofreq_call/{sample}_sorted.vcf.gz"
	output:
		single="lofreq_splitted/{sample}/{genome}.vcf",
		sorted="lofreq_splitted/{sample}/{genome}_sorted.vcf.gz",
		index="lofreq_splitted/{sample}/{genome}.vcf.tbi"
	shell: """
		tabix -h {input.vcf} {wildcards.genome} > {output.single}
		cat {output.single} | mawk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' | bgzip -c > {output.sorted}	 
		tabix -f -p vcf {output.sorted}
		cp {output.sorted}.tbi {output.index}  
		"""

rule split_bam:
	input:
		bam="mapped/{sample}.sorted.bam"
	output:
		single="mapped_splitted/{sample}/{genome}.bam"
	shell: """
		samtools view -bh {input.bam} > {output.single}
		samtools index {output.single}
		"""
rule annotate_vcf:
	input:
		vcf="lofreq_splitted/{sample}/{genome}_sorted.vcf.gz",
		bed="reference_genomes/{genome}.bed.gz"
	output:
		vcf2="lofreq_annotated/{sample}/{genome}.vcf"
	shell: """
		bcftools annotate -a {input.bed} -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') {input.vcf} > {output.vcf2}
		"""

rule snpEff:
	input:
		vcf="lofreq_splitted/{sample}/{genome}.vcf"
	params:
		ref="{genome}"
	output:
		vcf="snpEff/{sample}/{genome}.vcf",
		bed="snpEff/{sample}/{genome}.bed",
		csv="snpEff/{sample}/{genome}.csv",
		html="snpEff/{sample}/{genome}.html",
	log: 
		"logs/smpEff/{sample}/{genome}.log"
	shell: """				
		java -Xmx4g -jar ../tools/snpEff/snpEff.jar -no-downstream -no-upstream {params.ref} {input.vcf} -o bed > {output.bed}
		java -Xmx4g -jar ../tools/snpEff/snpEff.jar -no-downstream -no-upstream {params.ref} {input.vcf} -csvStats {output.csv} > {output.vcf} 2> {log}
		java -Xmx4g -jar ../tools/snpEff/snpEff.jar -no-downstream -no-upstream {params.ref} {input.vcf} -stats {output.html} > {output.vcf} 2> {log}
		""" 

rule igv_report:
	input:
		fasta="reference_genomes/{genome}.fasta",
				#vcf="snpEff/{sample}/{genome}.vcf",
		vcf="lofreq_splitted/{sample}/{genome}.vcf",
		bam="mapped_splitted/{sample}/{genome}.bam",
		bed="reference_genomes/{genome}.bed.gz"
	output:
		path="igv_report/{sample}/{genome}.html"
	shell: 
		"create_report {input.vcf} {input.fasta} --info-columns ANN --tracks {input.bam} {input.bed} --output {output.path}"

