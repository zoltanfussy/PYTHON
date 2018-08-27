#!/usr/bin/python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#faster than kika`s version, but relies on BioPython`s translate function, therefore a genetic code must be specified otherwise 
#ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt

infile = SeqIO.parse('/Volumes/zoliq data/OwnCloud/genomes/euglena gracilis/eghampl.fasta', 'fasta')
output = open('/Volumes/zoliq data/OwnCloud/genomes/euglena gracilis/eghampl.6frames.fasta', 'w')

i = 0
for sequence in infile:
    i += 1
    name = sequence.name
    seq = sequence.seq.upper()
    ambiguous = False
    for nucleotide in seq:
        if nucleotide not in 'ATCGN':
            ambiguous = True
    if ambiguous == True:
        with open('error.fasta', 'a') as error:
            error.write('>{}\n{}\n'.format(name, seq))
    else:    
        output.write('>{}_1\n{}\n'.format(name, str(seq.translate())))
        output.write('>{}_2\n{}\n'.format(name, str(seq[1:].translate())))
        output.write('>{}_3\n{}\n'.format(name, str(seq[2:].translate())))
        output.write('>{}_4\n{}\n'.format(name, str(seq.reverse_complement().translate())))
        output.write('>{}_5\n{}\n'.format(name, str(seq.reverse_complement()[1:].translate())))
        output.write('>{}_6\n{}\n'.format(name, str(seq.reverse_complement()[2:].translate())))
    if i%100 == 0:
        print(i)
output.close()