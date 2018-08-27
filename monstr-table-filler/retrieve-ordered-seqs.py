from Bio import SeqIO

to_retrieve = []
with open("temp.tsv") as infile:
	for line in infile:
		to_retrieve.append(line.strip())

infile_el = SeqIO.parse("el_merged.fasta", 'fasta')
my_nt_records = {}
for seq in infile_el:
	handle = seq.name
	if handle in to_retrieve:
		my_nt_records[handle] = seq
"""
infile_el = SeqIO.parse("../EGALL.nt.fasta", 'fasta')
for seq in infile_el:
	handle = seq.name.split(".")[0]
	if handle in to_retrieve:
		my_nt_records[handle] = seq
"""
ordered_records = []
for item in to_retrieve:
	ordered_records.append(my_nt_records[item])

SeqIO.write(ordered_records, "temp2-tobowtie.fasta", "fasta") #prints just the last seq, maybe a generator would help