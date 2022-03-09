from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import os
import argparse

print("This script uses a query fasta dataset, BLASTS it against a given database\nand analyses the output so that presequences (compared to best hits)\nare extracted and subjected to signalP and TMHMM.\n") 
print("Those sequences that have a SP or a TMD in their presequence are stored.\nSequences with no hit to the reference database will be shown via standard output.\n\nPlease make sure tmhmm and signalp are running in your system")

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "progs/PYTHON/N_extension_finder_Davca/"
os.chdir(wd)

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-db', '--database', help='DB of the huntee', required=True)
parser.add_argument('-q', '--query', help='fasta of genes to get', required=True)
parser.add_argument('-e', '--evalue', help='e-value threshold', default='1e-4')
parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='5')

args = parser.parse_args()

#makeblastdb -in input.fasta -out name_DB -dbtype prot -parse_seqids
if os.path.isfile("%s.phr" % (args.database)) == False:
    os.system('makeblastdb -in %s -out %s -dbtype prot -parse_seqids' % (args.database, args.database))
else:
    print("skipping makeblastdb... database exists")

print("db made")

#blastp -query 3plusTMD.fasta -db database.fasta -out blast.out -outfmt 5 -max_target_seqs 5
if os.path.isfile("blast.out") == False:
    os.system('blastp -query %s -db %s -out blast.out -outfmt 5 -evalue %s -max_target_seqs %s' % (args.query, args.database, args.evalue, args.num_seqs))

print ("blast done")

blast_out = open('blast.out')
blast_records = NCBIXML.parse(blast_out)

sequence_dict = {}
for sequence in SeqIO.parse(args.query, 'fasta'): # 3plusTMD.fasta
    sequence_dict[sequence.name] = sequence.seq



pandas_dict = {}
for record in blast_records:
    try:
        n = 0
        values_dict = {'start': [], 'hit_start': [], 'e_value': [], 'hit_id': []}
        for alignment in record.alignments:
            n += 1
            n_start = 1000000
            hit_id = alignment.hit_id
            for hsp in alignment.hsps:
                hsp_start = hsp.query_start
                if int(hsp_start) < n_start:
                    n_start = int(hsp_start)
                    hit_start = hsp.sbjct_start
                    e_value = hsp.expect
            values_dict['start'].append(n_start)
            values_dict['hit_start'].append(hit_start)
            values_dict['e_value'].append(e_value)
            values_dict['hit_id'].append(hit_id)
        pandas_dict[record.query] = pd.DataFrame(values_dict)
    except ValueError: 
        print(record.query)

with open('first_level', 'w') as res:
    for name, df in pandas_dict.items():
        try:
            best_blast_hit = df.loc[df['e_value'].idxmin()]['hit_id']
            med = int(df['start'].median())
            maxx = int(df['start'].max())
            sequence_med =  str(sequence_dict[name][:med])
            sequence_max =  str(sequence_dict[name][:maxx])
            if len(sequence) > 14:
                res.write('>{}_median {}\n{}\n'.format(name, best_blast_hit, sequence_med))
                res.write('>{}_max {}\n{}\n'.format(name, best_blast_hit, sequence_max))
        except ValueError: 
            print(name)

print("Waiting for predictions...")
os.system('signalp -t euk -f long first_level > first_level.signalp')
print("signalP successful")

signalp_dict = {}
for line in open('first_level.signalp'):
    if line.startswith('Name='):
        name = line.split()[0][5:]
        SP = line.split()[1]
        signalp_dict[name] = SP

#os.system('tmhmm first_level > first_level.tmhmm')
print("tmhmm successful")

tmhmm_dict = {}
for line in open('first_level.tmhmm'):
    if 'Number of predicted TMHs' in line:
        name = line.split()[1]
        TMHs = line.split()[6]
        tmhmm_dict[name] = TMHs

sequences = SeqIO.parse('first_level', 'fasta')

filtered = open('filtered_result.fasta', 'w')

with open ('final_result.fasta', 'w') as result:
    for sequence in sequences:
        result.write('>{} {} TMHs: {} BB_hit: {}\n{}\n'.format(sequence.name, 
                     signalp_dict[sequence.name], tmhmm_dict[sequence.name],
                     sequence.description.split()[1], sequence.seq))
        if 'max' in sequence.name:
            fullseqname = '_'.join(sequence.name.split('_')[:-1])
            if int(tmhmm_dict[sequence.name]) >= 1 or signalp_dict[sequence.name] == 'YES':
                filtered.write('>{} {} TMHs: {} BB_hit: {}\n{}\n'.format(fullseqname, 
                         signalp_dict[sequence.name], tmhmm_dict[sequence.name],
                         sequence.description.split()[1], sequence_dict[fullseqname]))

filtered.close()
print("Filtered sequences written, blasting filtered results ...")

os.system('blastp -query filtered_result.fasta -db %s -out filtered_blastp.out -outfmt 6 -evalue %s -max_target_seqs %s' % (args.database, args.evalue, args.num_seqs))

infile = open('filtered_blastp.out')
lines = infile.readlines()
infile.close()

#create a dictionary of lists
blast_d = {}
for line in lines:
    key = line.split('\t')[0]
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

for seq in seqs:
    seq_d[seq.split()[0]] = ''.join(seq.split('\n')[1:])

for key in blast_d:
    counter = 0
    out = open('%s.fas' % (key), 'w')
    for item in blast_d[key]:
        counter = counter + 1
        out.write('>%s\n%s\n' % (item, seq_d[item]))
    out.close()

print("Best hits to filtered sequences written, quitting.")

quit()













