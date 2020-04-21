'''
snakefile
by ,Philipp Muench

Maps NGS data to OligoMM reference genomes and performs variant calling
using the pooled samples
-----------------------------------------------------------------------

Dependencies:

samtools
	http://www.htslib.org/download/

fastp
	https://github.com/OpenGene/fastp

kneaddata
	https://bitbucket.org/biobakery/kneaddata

Trimmomatic
	http://www.usadellab.org/cms/?page=trimmomatic

bowtie2
	http://bowtie-bio.sourceforge.net/bowtie2/index.shtml	

Rscript
	https://cran.r-project.org/mirrors.html

sankeyD3 (R package)
	https://github.com/fbreitwieser/sankeyD3

networkD3 (R package)
	https://cran.r-project.org/web/packages/networkD3/index.html

bowtie2
	http://bowtie-bio.sourceforge.net/bowtie2/	

Usage:

	snakemake --jobs 20
'''

# Globals ---------------------------------------------------------------------

# name of sample to be processed
SAMPLES = ["1423_S34"]

FILES_CLAUDIA, = glob_wildcards("/net/metagenomics/projects/pmuench_oligomm_ab/ftp_reads/genome.gbf.de/HZI_BIFO/bifouser/19-0123/{files}_R1_001.fastq.gz")

# path to reads where .fastq.gz files are located
READS_PATH="/net/metagenomics/projects/pmuench_oligomm_ab/ftp_reads/genome.gbf.de/HZI_BIFO/bifouser/19-0123"

# dependencies
BOWTIE2_DIR="/net/metagenomics/projects/pmuench_oligomm_ab/tools/bowtie2-2.3.5.1-linux-x86_64"
PICARD_PATH="/net/metagenomics/projects/pmuench_oligomm_ab/tools/picard.jar"
TRIMMOMATIC_DIR="/net/metagenomics/projects/pmuench_oligomm_ab/tools/Trimmomatic-0.39"

# DBs
KRAKEN_DB="/net/sgi/oligomm_ab/tools/kraken2/kraken2-2.0.7-beta/db/minikraken2_v2_8GB_201904_UPDATE"
MEGARES_REFERENCE="/net/sgi/oligomm_ab/tools/megares/megares_drugs_database_v2.00.fasta" # ABR database reference
MEGARES_ANNOT="/net/sgi/oligomm_ab/tools/megares/megares_drugs_annotations_v2.00.csv" # ABR annotation

# bwa-indexed joined oligoMM reference file
REF="/net/sgi/oligomm_ab/OligoMM-report/databases/omm/joined_reference.fasta"
REF_ECOLI="/net/sgi/oligomm_ab/OligoMM-report/omm_plus_ecoli/joined_reference.fasta"

# Target rule -------------------------------------------------------------------

rule all:
	input:
		expand("docs/reports/fastp/{sample}.html", sample=SAMPLES),
#		expand("mapped/{sample}i.bam", sample=SAMPLES),
#		expand("vcf_varscan/{sample}.vcf", sample=FILES_CLAUDIA),
#		expand("lofreq/{sample}.vcf", sample=FILES_CLAUDIA),
#		expand("lofreq/{sample}.sorted.vcf.gz", sample=FILES_CLAUDIA),
#		expand("coverage/{sample}.coverage.txt.gz", sample=FILES_CLAUDIA),
#		expand("mapped_spacers/{sample}_1.bam", sample=SAMPLES),
#		expand("crispr_tables/{sample}.txt", sample=SAMPLES)
#		expand("megares/genes/{sample}.genes.txt", sample=FILES_PHILIPP),
#		expand("coverage_ind_plots/{sample}.{genome}.pdf", sample=FILES_PHILIPP, genome=GENOMES)

# Rules -------------------------------------------------------------------------

include: "rules/preprocess.smk" # copy files
#include: "../rules/varscan.smk" # variant calling 
#include: "../rules/bwa.smk" # mapping per sample
#include: "../rules/lofreq.smk" # SNP calling per genome and sample
#include: "../rules/coverage.smk" # coverage per sample per genome

