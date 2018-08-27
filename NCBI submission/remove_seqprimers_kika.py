#!/usr/bin/env python3
import os
from Bio import SeqIO
import re
#from collections import OrderedDict

with open('report.txt') as f:
	contamination = f.read().split("\n")
transcriptome = SeqIO.parse('el_merged_1line.fasta', 'fasta')
result = open('EL_merged_withoutNs_longer200_without_primers.fsa', 'w')

contigpattern = re.compile(r'Sequence-id: (\w+)')
primerpattern = re.compile(r'Interval: ([\.\w]+),')
primers = {}
for line in contamination:
	#print(line)
	if len(line) > 0:
		contig = contigpattern.search(line).group(1)
		primerhit = primerpattern.search(line).group(1)
		print(primerhit)
		start = int(primerhit.split("..")[0]) #POZOR neresi pripady kde jsou adaptery na obou koncich
		stop = int(primerhit.split("..")[1])
		primers[contig] = (start, stop)

for contig in transcriptome:
	key = contig.name
	if str(contig.seq).startswith('N'):
		pos = str(contig.seq).rfind('N')
		sequence = contig.seq[pos+1:]
	elif str(contig.seq).endswith('N'):
		pos = str(contig.seq).find('N')
		sequence = contig.seq[:pos]
	else:
		sequence = contig.seq

	if key in primers:
		if primers[key][0] == 1:
			sequence = sequence[primers[key][1]:]
		else:
			sequence = sequence[:primers[key][0]]
	if len(sequence) > 199:
		result.write('>{}\n{}\n'.format(key, sequence))


result.close()