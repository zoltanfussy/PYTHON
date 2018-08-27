#!/usr/bin/python
from Bio import SeqIO

infile = SeqIO.parse('/home/kika/blastocrithidia/datasets/Lpyr_CDs.fa', 'fasta')
taa = 0
tag = 0
tga = 0
other = 0

for sequence in infile:
	seq = sequence.seq[-3:]
	name = sequence.name
	if seq == 'taa':
		taa = taa+1
	elif seq == 'tag':
		tag = tag+1
	elif seq == 'tga':
		tga = tga+1
	else:
		other = other+1

stops = '{}: {}\n{}: {}\n{}: {}\n{}: {}'.format('taa',taa,'tag',tag,'tga',tga,'other',other)

with open('/home/kika/blastocrithidia/datasets/Lpyr_stops_number.txt', 'w') as result:
	result.write(stops)