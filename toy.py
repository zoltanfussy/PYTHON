"""with open('filtered.fasta') as f:
	line = f.read()
	seqs = line.split('>')[1:]

with open('translated.fasta', 'w') as out:
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
"""
my_big_list = [2, 6, 23, 54, 201]

if any(x < 3 for x in my_big_list):
	print("yeah")

