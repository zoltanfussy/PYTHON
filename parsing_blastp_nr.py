import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_protein
import pandas as pd
import argparse
import re

#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-d', '--directory', help='Set working directory', default='.')

args = parser.parse_args()

##################################
#### Create working directory ####
##################################

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
if args.directory == ".":
	print("changing to default directory")
	defdir = "Jankoviny/Tick_transcriptome/"
	wd = home + defdir
	os.chdir(wd)
else:
	os.chdir(args.directory)

#FILENAMES:
BLASTOUT = 'Trinity_all_trimmed_stages-AA_blastp2NR_db'
QUERYFILE = 'Trinity_all_trimmed_stages-AA.fasta'
FILTHITS = 'NR_filtered_hits.txt'
FILTFASTA = 'NR_goodproteins.fasta'
LOG = 'NR_seqextract-log.txt'
BLASTFETCHFILE = "blastfetch.txt"


#with open('example-BLAST.txt') as infile:
with open(BLASTOUT) as infile:
	lines = infile.read()
entries = lines.split('# BLASTP 2.7.1+\n')[1:]
print("{} entries found".format(len(entries)))
pandas_dict = {}

evaluethresh = 0.0001
goodqueries = set()
badqueries = set()
goodhits = set()
outfile = open(FILTHITS, 'w')
fetchfile = open(BLASTFETCHFILE, 'w')
c = 0
print("processing BLASTP output...")
for entry in entries:
	c += 1
	if c % 10000 == 0:
		print("{} queries processed".format(c))
	if len(entry) == 0:
		pass
	elif "# 0 hits found" in entry:
		pass
	else:
		entry = entry.split('\n') 
		query = entry[0].replace('# Query: ', '')
		#if query.startswith("c9811_g1_i1"):
		#	print("processing test contig analysis", query)
		queryname = query.split()[0]
		goodqueries.add(queryname)
		querylen = 5000 #int(re.search(r"len=([0-9]+)", query).group(1))
		hitsnumber = entry[3].replace('# ', '').replace(' hits found', '') 
		outfile.write("{}\t{}\nhit name\t% identity\talignment length\te-value\tq.start\tq.end\tq.frame\n".format(query, hitsnumber)) 
		counter = 0
		values_dict = {'hit_id': [], 'hit_start': [], 'e_value': [], 'absqstart': [], 'homologylen': [], 'querylen': []}
#		values_dict = {'hit_id': [], 'hit_start': [], 'e_value': [], 'align_len' : [], 'perc_ident': [], 'querylen': [], 'qstart': [], 'absqstart': [], 'qend': [], 'qframe': []}
		for hit in entry[4:]:
			counter += 1
			if len(hit.split()) > 0 and counter <= 10: #essentially top x hits
				hit_id = hit.split('\t')[1]
				hit_start = int(hit.split('\t')[8])
				e_value = float(hit.split('\t')[10])
				#print(e_value)
				if e_value > evaluethresh:
					badqueries.add(queryname)
				align_len = hit.split('\t')[3]
				perc_ident = hit.split('\t')[2]
				qstart = int(hit.split('\t')[6])
				qend = int(hit.split('\t')[7])
				if qstart < qend:
					absqstart = qstart
					homologylen = qend - qstart
				elif qstart > qend:
					print("qstart larger than qend")
					homologylen = qstart - qend
					absqstart = querylen - qstart + 1
				outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(hit_id, perc_ident, align_len, e_value, qstart, qend))
				if hit_id not in goodhits:
					goodhits.add(hit_id)
					fetchfile.write("{}\n".format(hit_id))
				values_dict['hit_id'].append(hit_id)
				values_dict['hit_start'].append(hit_start)
				values_dict['e_value'].append(e_value)
#				values_dict['align_len'].append(align_len)
#				values_dict['perc_ident'].append(perc_ident)
#				values_dict['qstart'].append(qstart)
				values_dict['absqstart'].append(absqstart)
				values_dict['homologylen'].append(homologylen)
				values_dict['querylen'].append(querylen)
#				values_dict['qend'].append(qend)
		pandas_dict[queryname] = pd.DataFrame(values_dict)
		outfile.write("==================\n")
print("not passing e-value threshold {}: {} {}".format(evaluethresh, len(badqueries), badqueries))

seq_count = 0
sequence_dict = {}
#testseq = open("testseq.fasta", "w")
for sequence in SeqIO.parse(QUERYFILE, 'fasta'):
#for sequence in SeqIO.parse('testseq.fasta', 'fasta'):
	if sequence.name in goodqueries:
		seq_count += 1
		sequence_dict[sequence.name] = sequence.description, sequence.seq
#		testseq.write(">{}\n{}\n".format(sequence.description, sequence.seq))
print("Sequences with positive hits extracted from dataset: {}".format(seq_count))

print("Now to adjusting positive query ORFs...")
with open(LOG, 'w') as logfile, open(FILTFASTA, 'w') as proteins:
	c = 0
	for sequence in sequence_dict:
		c += 1
		if c % 2500 == 0:
			print("{} sequences adjusted".format(c))
		nonsensemutation = False
		# checking if query is more than 10 aa shorter than best hit, flagging accordingly
		i = pandas_dict[sequence]['hit_start'][0]
		j = pandas_dict[sequence]['absqstart'][0]
		k = pandas_dict[sequence]['homologylen'][0]
		fromto = [j,k]
		if i > j + 10:
			orfN = "partial"
		else:
			orfN = "full"
		logfile.write("Query: {}\nhit_start {}, absolute_qstart {}, N terminus {}\n".format(sequence, i*3, j, orfN))

		# finding encoded protein, first upstream
		stops = ["*"]
		aalist = []
		downstream = sequence_dict[sequence][1]
		aacount = 0
		while aacount < j:
			aa = downstream[:1]
			aacount += 1
			downstream = downstream[1:]
			aalist.append(aa)
			if aa in stops:
				aa = ''
				aalist = []
		if "M" in aalist:
			start = aalist.index("M")
			protein_seq = sum(aalist[start:], Seq("", generic_protein))
			logfile.write("Extracted {} aa with start upstream of hit homology ({}).\n".format(len(protein_seq), protein_seq))
			newfrom = fromto[0] - len(protein_seq) + 1

		else:
			protein_seq = sum(aalist, Seq("", generic_protein))
			logfile.write("Extracted {} aa, NO START upstream of hit homology ({}).\n".format(len(protein_seq), protein_seq))
			newfrom = fromto[0] - len(protein_seq) + 1

		# finding encoded protein, now entering conserved region with the best hit
		minimumquerylen = len(protein_seq) + k
		while len(protein_seq) < minimumquerylen and len(downstream) > 0:
			aa = downstream[:1]
			if aa in stops:
				protein_seq = protein_seq + "X"
				nonsensemutation = True
			else:
				protein_seq = protein_seq + aa
			downstream = downstream[1:]

		# finding encoded protein, downstream of the homology region
		while aa not in stops and len(downstream) > 0:
			aa = downstream[:1]
			protein_seq = protein_seq + aa
			downstream = downstream[1:]

		#print(protein_seq[-3:])
		if aa in stops:
			orfC = "full"
		else:
			orfC = "partial"
		logfile.write("extracted {} aa total. extracted protein: {}\n======\n".format(len(protein_seq), protein_seq))
		fromto = [newfrom, newfrom + len(protein_seq)]
		#print("{} processed".format(sequence))

		if nonsensemutation:
			flag = "in-frame non-sense codon"
		elif orfC == "full" and orfN == "full":
			flag = "complete"
		elif orfC == "full" and orfN == "partial":
			flag = "N-partial"
		elif orfC == "partial" and orfN == "full":
			flag = "C-partial"
		else:
			flag = "fragment"
		proteins.write(">{} FLAG={} RANGE={}-{}\n{}\n".format(sequence_dict[sequence][0], flag, fromto[0], fromto[1] - 1, protein_seq))

print("Adjusted query proteins written to files.")
print("Please run the following to retrieve hits from database:")
print("fetch_from_local_BLASTDB.py")
print('grep ">" blastfetch.fasta > blastheaders.txt')

