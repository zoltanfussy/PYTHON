from Bio import SeqIO

infile = SeqIO.parse("bico_6frame_translated.fasta", "fasta")

with open("bico_6frame_translated-X.fasta", "w") as out:
	for seq in infile:
		sequence = str(seq.seq).replace("*", "X")
		out.write(">{}\n{}\n".format(seq.name, sequence))
