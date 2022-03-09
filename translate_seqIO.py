#!/usr/bin/python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#faster than kika`s version, but relies on BioPython`s translate function, therefore a genetic code must be specified otherwise 
#ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt

infile = SeqIO.parse('OM-RGC_viruses.fasta', 'fasta')

i = 0
with open('error.fasta', 'a') as error, open('OM-RGC_viruses.faa', 'w') as output:
	for sequence in infile:
		i += 1
		name = sequence.name
		seq = sequence.seq.upper()
		ambiguous = False
		for nucleotide in seq:
			if nucleotide not in 'ATCGN':
				ambiguous = True
				break
		if ambiguous == True:
				error.write('>{}\n{}\n'.format(name, seq))
		else:	
			output.write('>{}_1\n{}\n'.format(name, str(seq.translate())))
		if i%1000 == 0:
			print(i)
