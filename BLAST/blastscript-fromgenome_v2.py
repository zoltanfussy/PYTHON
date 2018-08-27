import os
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#this script uses tblastn to find hits in genomic data, then extracts the region, elongates the hit to nearest start/stop and outwrites protein translations
#the easy way for doing this is to use a all-reading-frame translated genomic database

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-db', '--database', help='DB of the huntee', required=True)
parser.add_argument('-q', '--query', help='fasta of genes to get', required=True)
parser.add_argument('-e', '--evalue', help='e-value threshold', default='1e-4')
parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='10')

args = parser.parse_args()

#check if database present, best if you check all the files, or just the last file that is made - database creation might have failed last time someone tried
#makeblastdb -in input.fasta -out name_DB -dbtype prot -parse_seqids
if os.path.isfile("%s.nhr" % (args.database)) == False:
	os.system('makeblastdb -in %s -out %s -dbtype nucl -parse_seqids' % (args.database, args.database))
else:
	print("skipping makeblastdb... database exists")

print('db made')

#blastp -query query.fas -db db -out name.out -outfmt 6 -evalue %s -max_target_seqs
"""
 -outfmt <String>
   alignment view options:
     0 = pairwise,
     1 = query-anchored showing identities,
     2 = query-anchored no identities,
     3 = flat query-anchored, show identities,
     4 = flat query-anchored, no identities,
     5 = XML Blast output,
     6 = tabular,
     7 = tabular with comment lines,
     8 = Text ASN.1,
     9 = Binary ASN.1,
    10 = Comma-separated values,
    11 = BLAST archive format (ASN.1) 
"""
os.system('tblastn -query {} -db {} -out temp.tblastn -outfmt 6 -evalue {} -max_target_seqs {}'.format(args.query, args.database, args.evalue, args.num_seqs))

print ('blast done')

infile = open('temp.tblastn')
lines = infile.readlines()
infile.close()

#create a dictionary of gene sets - contain not only hit name, but also start and end of hit - to enable extraction from genomic scaffolds
blast_d = {}
for line in lines:
	line = line.split('\t')
	key = line[0].split('_')[1]
	hit = line[1]
	start = line[8]
	end = line[9]
	if key in blast_d:
		blast_d[key].append((hit, start, end))
	else:
		blast_d[key]=[(hit, start, end)]

infile = open(args.database)
lines = infile.read()
infile.close()
seqs = lines.split('>')[1:]

seq_d = {}

for seq in seqs:
	seq_d[seq.split()[0]] = ''.join(seq.split('\n')[1:])

for key in blast_d:
	counter = 0
	out = open('fromgenome%s.fas' % (key), 'w')
	for item in blast_d[key]:
		scaffold = item[0]
		start = int(item[1])
		end = int(item[2])
		if start < end:
			contig_seq = seq_d[scaffold]
			#blast indexes sequences differently than python, so don`t ask why end is not -1, but this works as it should
			nucl_seq = contig_seq[start-1:end]
#extend ORF to the first Met upstream of hit, make sure nucl_seq[:3] is actually three nts
#			print (nucl_seq[:3])
#finding starting M
			stop_codons = ['TAG', 'TAA', 'TGA']
			if nucl_seq[:3] == 'ATG':
				pass
			else:
				triplet = 'XXX'
				upstream = contig_seq[:start-1]
				while triplet != 'ATG' and triplet not in stop_codons and len(triplet) == 3:
					triplet = upstream[-3:]
					nucl_seq = triplet + nucl_seq
					upstream = upstream[:-3]
			print ('found upstream M for ' + item)
			

			if nucl_seq[-3:] in stop_codons:
				pass
			else:
				triplet = 'XXX'
				downstream = contig_seq[end:]
				while triplet not in stop_codons and len(triplet) == 3:
					triplet = downstream[:3]
					nucl_seq = nucl_seq + triplet
					downstream = downstream[3:]
			print ('found downstream stop')
			coding_dna = Seq(nucl_seq, IUPAC.ambiguous_dna)
			protein_seq = coding_dna.translate()
			print (protein_seq)
			
		elif start > end:
			contig_seq = seq_d[scaffold]
			nucl_seq = contig_seq[end-1:start]
			coding_dna = Seq(nucl_seq, IUPAC.ambiguous_dna)
			protein_seq = coding_dna.reverse_complement().translate()
		else:
			pass
		counter = counter + 1
		out.write('>Giardia_%s@%s\n%s\n' % (counter, scaffold, protein_seq))
	out.close()

