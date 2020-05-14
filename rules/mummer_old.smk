# nucmer generates nucleotide alignments between two mutli-FASTA input files. The out.delta output file lists the distance between insertions and deletions that produce maximal scoring alignments between each sequence. 
rule nucmer:
	input: 
		fasta="databases/omm/ind_fasta/{ref}.fasta",
		reference="databases/omm/joined_reference.fasta"
	output: 
		delta="databases/omm/nucmer_delta/{ref}.delta"
	threads: 8
	log:
		"logs/nucmer/{ref}.txt"
	shell:
		"nucmer {input.reference} {input.fasta}  --delta {output.delta} --threads {threads}"

# show-coords displays a summary of the coordinates, percent identity, etc. of the alignment regions.
rule showcoords:
	input: 
		delta="databases/omm/nucmer_delta/{ref}.delta"
	output:
		tsv="databases/omm/duplication_coords/{ref}.tsv"
	shell: 
		"show-coords {input.delta} -o -c -H -l > {output.tsv}"
		
# Outputs a list of structural differences for each sequence in the reference and query, sorted by position. 
rule showdiff:
	input: 
		delta="databases/omm/nucmer_delta/{ref}.delta"
	output:
		tsv="databases/omm/duplication_diff/{ref}.tsv"
	shell: 
		"show-diff {input.delta} -H > {output.tsv}"
