import os
import sys
import argparse


database = sys.argv[1]
query = sys.argv[2]

print("if you want some e-value, use ...... -e 1e-4 .....")
print("if you want more than 10 hits, use ...... -n 50 .....")
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-e', '--evalue', help='e-value threshold', default='1e-4')
parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='10')

args = parser.parse_args()

#makeblastdb -in input.fasta -out name_DB -dbtype prot -parse_seqids


if exists("/home/kika/programs/blast-2.5.0+/bin/{}.nhr" .format(database)) == False:
	print("database is going to be made")
	os.system('makeblastdb -in {} -out {} -dbtype nucl -parse_seqids' .format(database, database))
	print('db made')
else:
	print("skipping makeblastdb... database exists")

#blastp -query query.fas -db db -out name.out -outfmt 6 -evalue %s -max_target_seqs
#-outfmt 6 - table format, tab separated %s
os.system('tblastn -query {} -db {} -out temp.blastp -outfmt 5 -evalue {} -max_target_seqs {}' .format(query, database, args.evalue, args.num_seqs))

print('blast done')

infile = open('temp.blastp')
lines = infile.readlines()
infile.close()

