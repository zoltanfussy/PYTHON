from Bio import SeqIO
import re

infile = SeqIO.parse("uniprot_sprot", "fasta")
output = open('sprot_parsed.fasta', 'w')

#multiplications = {}
seq_dict = {}
unwanted = ['/', '(', ')', ',']
for sequence in infile:
	#multiplications[sequence.seq].append(sequence.name)
	if sequence.seq not in seq_dict:
		seqname = sequence.name
		if "tr|" in seqname:
			seqname = seqname.replace('tr|', '')
		if "sp|" in seqname:
			seqname = seqname.replace('sp|', '')
		taxcode = seqname.split('_')[1]
		acc = seqname.split('|')[0]
		genecode = seqname.split('|')[1].split('_')[0]
		seqname = ('{}@{}_{}'.format(genecode, taxcode, acc))
		"""print("ID: ", sequence.id)
		print("NAME: ", sequence.name)
		print("LETT_ANOT: ", sequence.letter_annotation)
		print("ANNOTS: ", sequence.annotations)
		print("FEATURES: ", sequence.features)
		print("dbxref: ", sequence.dbxrefs)"""
		seqdescription = sequence.description
		for item in unwanted:
			while item in seqdescription:
				seqdescription = seqdescription.replace(item, '')
		annotation = seqdescription.split(' OS=')[0].split(' ')[1:]
		annotation = ' '.join(annotation)

		taxon = seqdescription.split('OS=')[1].split(' GN=')[0]
		if 'isolate' in taxon:
			taxon = taxon.replace('isolate ', 'isol.')
		#print(seqname, taxon, annotation)
		fullannotation = ('{} {} [{}]'.format(seqname, annotation, taxon))
		seq_dict[sequence.seq] = fullannotation

for key, value in seq_dict.items():
	output.write('>{}\n{}\n'.format(value, key))

output.close()
