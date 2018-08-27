infile1 = open('infile.txt')
line1 = infile1.read()
infile1.close()

out1 = open('sed-key-sprot.txt','w')

seqs = line1.split('\n')

for sequence in seqs:
	writeseq = "s/" + sequence + "/g\n"
	out1.write(writeseq)

out1.close()