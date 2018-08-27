from Bio import SeqIO
import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()

files = os.listdir('.')
files = [f for f in files if f.endswith(".pep.all.fa") and os.path.isfile(f)]

with open("ryby_prot.accession2taxid", "w") as result:
	for file in files:
		c = 0
		taxname = [file.split(".")[0].replace("_"," ")]
		print(taxname)
		taxid = ncbi.get_name_translator(taxname)[taxname[0]][0]
		for sequence in SeqIO.parse(file, "fasta"):
			c += 1
			result.write("{0}\t{0}\t{1}\t000000000\n".format(sequence.name, taxid))
			if c % 1000 == 0:
				print(c)
