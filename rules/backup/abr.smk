# align PE reads to antibiotic database reference
rule bwa_mem_res:
	input:
		reads=["trimmed/pe/{sample}.1.fastq.gz", "trimmed/pe/{sample}.2.fastq.gz"]
	output:
		"mapped_megares/{sample}.bam"
	log:
		"logs/mapped_megares/{sample}.log"
	params:
		index=MEGARES_REFERENCE,
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort="none",						 # Can be 'none', 'samtools' or 'picard'.
		sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
		sort_extra=""						 # Extra args for samtools/picard.
	wrapper:
		"0.36.0/bio/bwa/mem"

# convert bam to sam
rule bam_to_sam:
				input:
								bam="mapped_megares/{sample}.bam"
				output:
								sam="mapped_megares/{sample}.sam"
				shell:
								"samtools view -h -o {output.sam} {input.bam}"

# process sam with megares
rule sam_to_megares:
				input:
								sam="mapped_megares/{sample}.sam"
				output:
								genes="megares/genes/{sample}.genes.txt",
								groups="megares/groups/{sample}.groups.txt",
								mech="megares/mech/{sample}.mech.txt",
								cla="megares/class/{sample}.class.txt"
				params:
								ref=MEGARES_REFERENCE,
								annot=MEGARES_ANNOT
				log:
								"logs/megares/{sample}.log"
				shell:  """
								resistome -ref_fp {params.ref} \
																-annot_fp {params.annot} \
																-sam_fp {input.sam} \
																-gene_fp {output.genes} \
																-group_fp {output.groups} \
																-mech_fp {output.mech} \
																-class_fp {output.cla} \
																-t 80
								"""

