#function to read fasta and create a list of sequences
def read_fasta(filename):
	infile = open(filename)
	lines = infile.read()
	sequences = lines.split(">")[1:]
	infile.close()
	return sequences

#function to count GC in a single sequence
def GC_count(sequence):
	G_count = sequence.count('G') + sequence.count('g')
	C_count = sequence.count('C') + sequence.count('c')
	A_count = sequence.count('A') + sequence.count('a')
	T_count = sequence.count('T') + sequence.count('t')
	GC_perc = 100*(float(G_count + C_count) / (G_count + C_count + A_count + T_count))
	return GC_perc

#exercise function to reverse complement a sequence
def reverse_complement(sequence):
	rc_seq = ""
	# make dictionary of reverse complements
	revcomplements = {}
	revcomplements['A'] = 'T'
	revcomplements['a'] = 'T'
	revcomplements['T'] = 'A'
	revcomplements['t'] = 'A'
	revcomplements['C'] = 'G'
	revcomplements['c'] = 'G'
	revcomplements['G'] = 'C'
	revcomplements['g'] = 'C'

	for base in sequence:
#instead of writing upper and lower case bases, you can uppercase the variable
		base = base.upper()
		if base in revcomplements:
			rc_base = revcomplements[base]
		else:
			rc_base = 'N'
		rc_seq = rc_base + rc_seq
	return rc_seq	

seqs = read_fasta('sequence.fasta')
outfile = open('rc_sequence.fasta','w')
for seq in seqs:
	seqname = seq.split('\n')[0]
	full_sequence = ''.join(seq.split('\n')[1:])
	if len(full_sequence) <= 20000:
		rc_seq = reverse_complement(full_sequence)
		outfile.write('>'+ seqname + '\n' + rc_seq + '\n')
	else:
		pass
outfile.close()