
rule snpEff:
	input:
		vcf="lofreq_splitted/{sample}_{genome}.vcf"
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


