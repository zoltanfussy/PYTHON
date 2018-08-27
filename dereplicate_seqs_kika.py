#from collections import OrderedDict

from Bio import SeqIO

infile = SeqIO.parse('FNR_sequences_v2.fasta', 'fasta')
output1 = open('outfile_dedupl_FNR.fasta', 'w')
output2 = open('outfile_dupl_names.fasta', 'w')

multiplications = {}
#seq_dict = OrderedDict()
seq_dict = {}
for sequence in infile:
	if sequence.seq in multiplications:
		multiplications[sequence.seq].append(sequence.name)
	else:
		multiplications[sequence.seq] = [sequence.name]
	if sequence.seq not in seq_dict:
		#rename full header only with name (acc, till the first space)
		# seq_dict[sequence.seq] = sequence.name 
		#keep full header
		seq_dict[sequence.seq] = sequence.description

for key, value in seq_dict.items():
	output1.write('>{}\n{}\n'.format(value, key))

for key, value in multiplications.items():
	if len(value) > 1:
		output2.write('{}\n'.format(str(value)))

output1.close()
output2.close()