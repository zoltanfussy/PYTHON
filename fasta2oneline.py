from Bio import SeqIO

file = SeqIO.parse("el_merged.fasta", "fasta")

with open("el_merged_1line.fasta", "w") as result:
	for seq in file:
		result.write(">{}\n{}\n".format(seq.description, seq.seq))