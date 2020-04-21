# run crass on fastq

rule crass:
  input:
    fq="trimmed/pe/{sample}.1.fastq.gz"
  output:
    crass=directory("crass/{sample}_1")
  shell: """
	       /net/sgi/oligomm_ab/tools/crass/crass-0.3.12/bin/bin/crass -o {output.crass} {input.fq}
				 """
rule export_spacers:
  input:
    crispr="crass/{sample}_1/crass.crispr"
  output:
    txt="crass/{sample}_1/crispr.txt"
	shell: """
         grep -e "<spacer cov=" -e  "<group drseq=" {input.crispr} > {output.txt}
         """

rule export_spacers2:
  input:
    crispr="crass/{sample}_1/crass.crispr"
  output:
    txt="crispr_sequence/{sample}.txt"
  shell: """
         grep -e "<spacer cov=" -e  "<group drseq=" {input.crispr} > {output.txt}
         """

rule spacertable:
  input:
    crass="crass/{sample}_1/crass.crispr"
  output:
    csv="crispr_tables/{sample}.txt"
  shell: """
         python /net/sgi/oligomm_ab/tools/crass2table.py {input.crass} {output.csv} {input.crass}
         """

rule spacer_fasta:
  input:
    txt="crass/{sample}_1/crispr.txt"
  output:
    fasta="crass/{sample}_1/spacers.fasta"
  shell: """
	       grep "spacer cov" {input.txt} |  awk -F '"' '{{print $4}}'  | awk '{{print ">"NR"\\n"$0}}' > {output.fasta}
	"""

rule spacer_bowtie2:
  input:
    sample="crass/{sample}_1/spacers.fasta"
  output:
    bam="mapped_spacers/{sample}_1.bam"
  shell: "bowtie2 -x /net/sgi/oligomm_ab/databases/oligomm_full_bt2/ref -f {input.sample} --no-unal |  samtools view -S -b > {output.bam}"
