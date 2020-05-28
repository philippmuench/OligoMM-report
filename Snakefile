"""Based on the best-practices variant calling implementation LoFreq
which can be found here https://github.com/gis-rpd/pipelines and in the
lofreq github repo. Finishes with a bgzipped vcf file.

# Input: config-file with following fields:
- bool mark_short_splits: for bwa mem -M
- string bed: for bed-file limiting analysis to certain regions
- int optional 'maxdepth': for limit per-site coverage in analysis
- dict 'samples': sample names as keys and one fastq-pair each as value
- string reference: reference fasta file
- string outdir: where to save output

# Pre-installed programs:
- lofreq 2.1.2
- bwa (with mem support e.g. 0.7.12)
- samtools >= 1.3
- bamCoverage
- pilon

Notes:
- If missing, the workflow will try to index your reference with 
  samtools and bwa. This can lead to race conditions so is best
  done in advance.
"""

import os

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail;")

##########################################################################
# define output files		   		      	 		 #
##########################################################################

rule all:
    input:
        expand(os.path.join(config['outdir'],
	"{sample}/{sample}.vcf.gz"),
	sample=config['samples'])

# helpers
def get_file(s, elem = 0, trimmed = False):
    if trimmed:
        t = [sub + '.trimmed.fastq.gz' for sub in s]
        return t[elem]
    else:
        return s[elem]

##########################################################################
# define rules   		   		      	 		 #
##########################################################################

rule bwa_index:
    input:
        "{prefix}.{suffix}"
    output:
        "{prefix}.{suffix,(fasta|fa)}.pac",
        "{prefix}.{suffix,(fasta|fa)}.bwt",
        "{prefix}.{suffix,(fasta|fa)}.sa"
    log:
        "{prefix}.{suffix,(fasta|fa)}.index.log"
    shell:
        "bwa index {input} >& {log};"

rule samtools_faidx:
    input:
        "{prefix}.{suffix}"
    output:
        "{prefix}.{suffix,(fasta|fa)}.fai",
    log:
        "{prefix}.{suffix,(fasta|fa)}.index.log"
    shell:
        "samtools faidx {input} >& {log};"

rule fastp:
    input:
        fq1 = lambda wc: get_file(config['samples'][wc.sample], elem = 0),
        fq2 = lambda wc: get_file(config['samples'][wc.sample], elem = 1),
    output:
        html = "{prefix}/{sample}.html",
        o1 = temp("{prefix}/{sample}.R1.fastq.gz"),
        o2 = temp("{prefix}/{sample}.R2.fastq.gz")
    log:
        "{prefix}/{sample}.fastp.log"
    threads:
        8
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2}"
        " --out1 {output.o1} --out2 {output.o2} --html {output} >& {log};"

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai",
    log:
        "{prefix}.bam.bai.log"
    shell:
        "samtools index {input} >& {log};"
        
rule bam_coverage:
    input:
         bam="{prefix}.bam",
         bai="{prefix}.bam.bai"
    output:
        "{prefix}.bam.bw",
    log:
        "{prefix}.bam.bw.log"
    shell:
        "bamCoverage -b {input.bam} -o {output} >& {log};"

rule bwamem_align:
    input:
        reffa = lambda wc: config['references'][wc.sample], 
        bwaindex = lambda wc: config['references'][wc.sample] + ".bwt",
        fq1 = "{prefix}/{sample}.R1.fastq.gz",
        fq2 = "{prefix}/{sample}.R2.fastq.gz"
    output:
        bam = temp("{prefix}/{sample}.bwamem.bam")
    log:
        "{prefix}/{sample}.bwamem.bam.log"
    params:
        mark_short_splits = "-M" if config['mark_short_splits'] else "",
    message:
        "Aligning PE reads, fixing mate information and converting to sorted BAM"
    threads:
        8
    shell:
        "{{ bwa mem {params.mark_short_splits} -t {threads}"
        " {input.reffa} {input.fq1} {input.fq2} |"
        " samtools fixmate - - |"
	" samtools view -bq 1 |"
        " samtools sort -o {output.bam} -T {output.bam}.tmp -; }} >& {log}" 

rule bwamem_unfiltered:
    """ Generates unfiltered .bam file for comparison in NBG, not included as a
    output in rule_all
    """
    input:
        reffa = lambda wc: config['references'][wc.sample], 
        bwaindex = lambda wc: config['references'][wc.sample] + ".bwt",
        fq1 = "{prefix}/{sample}.R1.fastq.gz",
        fq2 = "{prefix}/{sample}.R2.fastq.gz"
    output:
        bam = "{prefix}/{sample}.unfiltered.bam"
    log:
        "{prefix}/{sample}.unfiltered.bam.log"
    params:
        mark_short_splits = "-M" if config['mark_short_splits'] else "",
    message:
        "Aligning PE reads, fixing mate information and converting to sorted BAM"
    threads:
        8
    shell:
        "{{ bwa mem {params.mark_short_splits} -t {threads}"
        " {input.reffa} {input.fq1} {input.fq2} |"
        " samtools sort -o {output.bam} -T {output.bam}.tmp -; }} >& {log}" 

rule lofreq_bam_processing:
    """Runs BAM through full LoFreq preprocessing pipeline,
    i.e. viterbi, alnqual, indelqual, followed by sort (required by
    viterbi).

    WARNING: running this on unsorted input files will be inefficient
    because of constant reloading of the reference
    """
    input:
        bam = "{prefix}/{sample}.bwamem.bam",
        reffa = lambda wc: config['references'][wc.sample],
        reffai = lambda wc: config['references'][wc.sample] + ".fai",
        refadditionalbam = "{prefix}/{sample}.bwamem.bam.bai"
    output:
        bam = '{prefix}/{sample}.lofreq.bam'
    log:
        '{prefix}/{sample}.lofreq.log'
    message:
        "Preprocessing BAMs with LoFreq"
    threads:
        1
    shell:
        "{{ lofreq viterbi -f {input.reffa} {input.bam} | "
        " lofreq alnqual -u - {input.reffa} | "
        " lofreq indelqual --dindel -f {input.reffa} - | "
        " samtools sort -o {output.bam} -T {output.bam}.tmp -; }} >& {log}"

rule pilon_refinement:
    """Runs pilon on reference seqeunce using bam
    """
    input:
        bam = "{prefix}/{sample}.bwamem.bam",
        bai = "{prefix}/{sample}.bwamem.bam.bai",
        reffa = lambda wc: config['references'][wc.sample],
        reffai = lambda wc: config['references'][wc.sample] + ".fai",
    output:
        pilon = "{prefix}/{sample}_pilon/",
        changes = "{prefix}/pilon/{sample}.changes",
    log:
        '{prefix}/{sample}_pilon/pilon.log'
    message:
        "run pilon on reference"
    threads:
        8
    params:
        output_dir = "out/{sample}/pilon",
        output_prefix = "{sample}",
        max_mem = "80G"
    shell:
        'pilon -Xmx{params.max_mem} --genome {input.reffa} --bam {input.bam}'
        ' --output {params.output_prefix} --outdir {params.output_dir}'
        ' --changes --fix all --threads {threads} 2> {log}'

rule lofreq_call:
    input:
        bam = "{prefix}/{sample}.lofreq.bam",
        bai = "{prefix}/{sample}.lofreq.bam.bai",
        coverage = "{prefix}/{sample}.lofreq.bam.bw",
        reffa = lambda wc: config['references'][wc.sample],
        refidx = lambda wc: config['references'][wc.sample] + ".fai"
    output:
        vcf = '{prefix}/{sample}.vcf.gz'
    log:
        '{prefix}/{sample}.vcf.log'
    message:
        "Calling variants with LoFreq"
    threads:
        8
    params:
        maxdepth = config.get('maxdepth', 10000),
        bed_arg = "-l {}".format(config['bed']) if config['bed'] else ""
    shell:
        "lofreq call-parallel --pp-threads {threads} --call-indels"
        " {params.bed_arg} -f {input.reffa} -o {output.vcf}"
        " -d {params.maxdepth} {input.bam} >& {log}"
