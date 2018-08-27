#!/usr/bin/python3
from Bio import SeqIO

infile = SeqIO.parse('./Trinity_Sterlet_mat,oos_WTA.fa', 'fasta')
output = open('./count.txt', 'w')

i = 0
for sequence in infile:
    i += 1

output.write(i)
output.close()