from Bio import SeqIO

seqset = set()
houby = SeqIO.parse("houby_derepped.fasta", "fasta")
for seq in houby:
	seqset.add(seq.seq)

files = [x for x in os.listdir(".") if x.endswith("ryby.fasta")]
with open("refseq_houby.fasta", "w") as derep:
	for file in files:
		seqdata = SeqIO.parse(file, "fasta")
		for seq in seqdata:
			if seq.seq not in seqset:
				#seqset.add(seq.seq)
				derep.write(">{}\n{}\n".format(newname, seq.seq))
