from Bio import SeqIO
import os

infile = SeqIO.parse("cleaved-vbra_FINALSETinclEXT.fasta", "fasta")
maxseqcount = 3300 #a little more than a quarter of total seq count

counter = 0
directory = 1
os.system("mkdir 1")
os.system("mkdir 2")
os.system("mkdir 3")
os.system("mkdir 4")

for sequence in infile:
	counter += 1
	with open("{}/{}.fasta".format(directory, sequence.name), "w") as outfile:
		outfile.write(">{}\n{}\n".format(sequence.description, sequence.seq))
	if counter == maxseqcount:
		directory += 1
		counter = 0
