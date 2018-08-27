import os
import argparse

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-db', '--database', help='DB of the huntee', required=True)
parser.add_argument('-q', '--query', help='fasta of genes to get', required=True)
parser.add_argument('-e', '--evalue', help='e-value threshold', default='1e-4')
parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='10')

args = parser.parse_args()

#makeblastdb -in input.fasta -out name_DB -dbtype prot -parse_seqids
if os.path.isfile("%s.nhr" % (args.database)) == False:
	os.system('makeblastdb -in %s -out %s -dbtype prot -parse_seqids' % (args.database, args.database))
else:
	print "skipping makeblastdb... database exists"

print 'db made'

#blastp -query query.fas -db db -out name.out -outfmt 6 -evalue %s -max_target_seqs
#-outfmt 6 - table format, tab separated %s
os.system('blastp -query %s -db %s -out temp.blastp -outfmt 6 -evalue %s -max_target_seqs %s' % (args.query, args.database, args.evalue, args.num_seqs))

print 'blast done'

infile = open('temp.blastp')
lines = infile.readlines()
infile.close()

#create a dictionary of lists
blast_d = {}
for line in lines:
	key = line.split('\t')[0]
	key = key.split('_')[1]
	hit = line.split('\t')[1]
	if key in blast_d:
		blast_d[key].append(hit)
	else:
		blast_d[key]=[hit]

infile = open(args.database)
lines = infile.read()
infile.close()
seqs = lines.split('>')[1:]

seq_d = {}

#one way to export fastas using two dictionaries
for seq in seqs:
	seq_d[seq.split()[0]] = ''.join(seq.split('\n')[1:])

for key in blast_d:
	counter = 0
	out = open('%s.fas' % (key), 'w')
	for item in blast_d[key]:
		counter = counter + 1
		out.write('>Giardia%s\n%s\n' % (counter, seq_d[item]))    # the sequence name is Giardia1, Giardia2, etc in the output file

	out.close()

#a second option is to create a list of hits and then run blastdbcmd to retrieve the sequences
#this second way there might be problems with sequence names