#import pandas as pd
import os
import argparse
import time
from Bio.Blast import NCBIXML
from Bio import SeqIO

#BACHA běží déle než verze 1, nicméně jsem nahradil allhits_list allhits_set, takže se to mohlo zrychlit
print("Start time:", time.ctime())
#BLAST setup:
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-db', '--database', help='DB of the huntee', required=True)
# parser.add_argument('-q', '--query', help='fasta of genes to get', required=True)
# parser.add_argument('-e', '--evalue', help='e-value threshold', default='1e-4')
# parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='5')

args = parser.parse_args()


# #makeblastdb -in input.fasta -out name_DB -dbtype prot -parse_seqids
# if os.path.isfile("%s.phr" % (args.database)) == False:
#     os.system('makeblastdb -in %s -out %s -dbtype prot -parse_seqids' % (args.database, args.database))
# else:
#     print("skipping makeblastdb... database exists")
# print("db made")

# #blastp -query 3plusTMD.fasta -db database.fasta -out blast.out -outfmt 5 -max_target_seqs 5
# if os.path.isfile("blast.out") == False:
#     os.system('blastp -query %s -db %s -out blast.out -outfmt 5 -evalue %s -max_target_seqs %s' % (args.query, args.database, args.evalue, args.num_seqs))
# print ("blast done")


# load query sequences
# sequence_dict = {}
# for sequence in SeqIO.parse(args.query, 'fasta'):
#     sequence_dict[sequence.name] = sequence.seq


#reading blast output
print("Reading BLAST output")
files = os.listdir()
hits_dict = {}
allhits = set()
for file in files:
    if "tblastn" in file:
        try:
            blast_out = open(file)
            blast_records = NCBIXML.parse(blast_out)
            for record in blast_records:
                record_id = record.query.split()[0]
                hits = []
                for alignment in record.alignments:
                    hit_id = alignment.hit_id
                    hits.append(hit_id)
                    allhits.add(hit_id)
                hits_dict[record_id] = hits
        except ValueError: 
            print(record.query)

# load hit sequences
print("Loading database sequences, might take a while...")
database_dict = {}
for sequence in SeqIO.parse(args.database, 'fasta'):
    if sequence.name in allhits:
        database_dict[sequence.name] = sequence.seq


print("Creating results directory")
if os.path.isdir("results") != True:
    os.mkdir("results")
os.chdir("results")

#writing of found hits
for key in hits_dict.keys():
    with open (key + '.fasta', 'w') as result:
        for hit in hits_dict[key]:
            result.write('>{}@{}\n{}\n'.format(key, hit, database_dict[hit]))

print("Found hits sequences written")

#writing statistics
with open('statistics.txt', 'w') as statistics:
    for key in hits_dict.keys():
        statistics.write("{}: {} reads found\n".format(key, len(hits_dict[key])))
print("# found reads written into statistics.txt")

print("Finish time:", ctime())

quit()
