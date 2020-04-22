rule igv_report_omm:
	input:
		fasta=config["REF"],
		vcf="processed/lofreq_bwa_sorted/{sample}.vcf.gz",
		# any number of additional optional tracks, see igv-reports manual
		tracks="processed/bwa_mapped_omm/{sample}.bam"
	output:
		"docs/igv/{sample}.html"
	params:
		extra=""  # optional params, see igv-reports manual
	log:
		"logs/igv/{sample}.log"
	wrapper:
		"0.51.3/bio/igv-reports"
