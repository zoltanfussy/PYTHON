from Bio import SeqIO

infile = SeqIO.parse("el_merged.fasta", "fasta")
n=0
m=0
for sequence in infile:
	if 'TTTTTCG' in sequence.seq[:35]:
		n += 1
	elif 'CGAAAAA' in sequence.seq[-35:]:
		m += 1
print(n)
print(m)
print(n+m)