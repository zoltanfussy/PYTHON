import os
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#this script uses tblastn to find hits in genomic data, then extracting

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
	print "skipping makeblastdb... database exists"

print 'db made'

#blastp -query query.fas -db db -out name.out -outfmt 6 -evalue %s -max_target_seqs
#-outfmt 6 - table format, tab separated %s
os.system('tblastn -query %s -db %s -out temp.tblastn -outfmt 6 -evalue %s -max_target_seqs %s' % (args.query, args.database, args.evalue, args.num_seqs))

print 'blast done'

infile = open('temp.tblastn')
lines = infile.readlines()
infile.close()

#create a dictionary of gene sets - contain not only hit name, but also start and end of hit - to enable extraction from genomic scaffolds
blast_d = {}
for line in lines:
	key = line.split('\t')[0]
	key = key.split('_')[1]
	hit = line.split('\t')[1]
	start = line.split('\t')[8]
	end = line.split('\t')[9]
	if key in blast_d:
		blast_d[key].append((hit, start, end))
	else:
		blast_d[key]=[(hit, start, end)]

infile = open(args.database)
lines = infile.read()
infile.close()
seqs = lines.split('>')[1:]

seq_d = {}

#one way to export fastas using two dictionaries
for seq in seqs:
	seq_d[seq.split()[0]] = ''.join(seq.split('\n')[1:])

for key in blast_d:
	print key, 'this is gene name'
	counter = 0
	out = open('fromgenome%s.fas' % (key), 'w')
	for item in blast_d[key]:
		scaffold = item[0]
		start = int(item[1])
		end = int(item[2])
		if start < end:
			contig_seq = seq_d[scaffold]
			#blast indexes sequences differently than python, so don`t ask why end is not -1, but this works
			nucl_seq = contig_seq[start-1:end]
			coding_dna = Seq(nucl_seq, IUPAC.ambiguous_dna)
			protein_seq = coding_dna.translate()
			print nucl_seq
			print protein_seq
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

#a second option is to create a list of hits and then run blastdbcmd to retrieve sequences
#this second way there might be problems with sequence names