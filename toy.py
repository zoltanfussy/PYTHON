infile = open('filtered.fasta')
line = infile.read()
infile.close()

seqs = line.split('>')[1:]

out = open('translated.fasta', 'w')

for seq in seqs:
	coding_dna = ''.join(seq.split('\n')[1:])
	if len(coding_dna) % 3 == 2:
		coding_dna = coding_dna[:-2]
	elif len(coding_dna) % 3 == 1:
		coding_dna = coding_dna[:-1]
	else:
		pass
	print len(coding_dna) % 3
	out.write(">%s\n%s\n" % (seq.split('\n')[0], coding_dna))
out.close()