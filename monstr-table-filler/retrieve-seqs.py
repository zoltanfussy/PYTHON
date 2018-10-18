from Bio import SeqIO

to_retrieve = set()
with open("nejakyTvujSeznam.tsv") as infile:
	for line in infile:
		to_retrieve.add(line.strip())

infile_el = SeqIO.parse("databaze.fasta", 'fasta')
with open("out.fasta", "w") as result:
	for seq in infile_el:
		if seq.name in to_retrieve:
			result.write(">{}\n{}\n".format(seq.name, seq.seq))
