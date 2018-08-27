from Bio import SeqIO,AlignIO
import os

files = os.listdir('.')
files = [f.replace(".fasta","") for f in files if f.endswith(".fasta")]
for file in files:
	with open("{}.fasta".format(file)) as infile, open("{}.phy".format(file), "w") as outfile:
		alignmentfile = AlignIO.read(infile, "fasta")
		AlignIO.write(alignmentfile, outfile, "phylip-relaxed")