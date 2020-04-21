# calculate coverage on deduplicated bam files
rule samtools_coverage:
		input: 
				 bam="mapped/{sample}.bam"
		output:
				 stat="coverage/{sample}.coverage.txt.gz"
		shell: "samtools depth -a {input.bam} | bgzip > {output.stat}"

#rule plot_coverage:
#		input:
#			stat="coverage_ind/{sample}.{genome}.coverage.txt.gz"
#		output:
#			pdf="coverage_ind_plots/{sample}.{genome}.pdf"
#		shell: "Rscript ../scripts/plotCoverageHist.R {input.stat} {output.pdf}"

