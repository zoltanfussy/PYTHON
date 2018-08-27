from Bio import SeqIO

seqset = set()
with open("CV&VB_heme_pathway.txt") as file:
	seqlist = file.readlines()
	for seq in seqlist:
		seqset.add(seq.split(" ")[0])

cvel = SeqIO.parse("chromera_all-elong_aa_mmetsp-ncbi.txt", "fasta")
vbra = SeqIO.parse("vitrella_all-elong_aa_mmetsp.txt", "fasta")

for seq in cvel:
	for item in seqset:
		if seq.name.startswith(item):
			print(">{}\n{}\n".format(seq.name, seq.seq))

for seq in vbra:
	for item in seqset:
		if seq.name.startswith(item):
			print(">{}\n{}\n".format(seq.name, seq.seq))