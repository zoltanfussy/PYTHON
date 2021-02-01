import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#os.system('source ~/.bash_profile')

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-in', '--infile', help='fasta to be processed', required=True)
parser.add_argument('-sl', '--SL', help='sequence of the spliced leader', default='TTTTTCG')

args = parser.parse_args()

filename = args.infile.split(".")[0]

monstr_dic = {}
###########################################################################
#        read in seq names to extract nt/aa models of interest:           #
###########################################################################
if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
defdir = "progs/PYTHON/monstr-table-filler/"
os.chdir(home + defdir)
 
names = []
with open(args.infile) as f:
	for l in f:
		line = l.strip()
		name = line
		name = name.split("_")[0]
		if "-reversed" in name:
			name = name.replace("-reversed", "")
		if name not in names:	
			names.append(name)
		monstr_dic[name] = 	["annotation", "abbrev", name, "3=reads", "4=length", "5=avg_cov", "6=SL", "7=TMHMM", "8=PrediSi", "9=PredSL", "no SP, TP N/D", "11=nt"]

#monstr_dic["title"] = ["annotation", "short", "contig", "no.reads", "contig len", "avg_cov", "SL", "TMHMM", "PrediSi", "PredSL", "TP_follow_SP", "nt seq", "aa seq"]
#monstr_dic["indices"] = [  "0",        "1",      "2",       "3",         "4",        "5",    "6",    "7",       "8",     "9",         "10",        "11",     "12"]


intranscripts = SeqIO.parse(home + "genomes/euglena longa/ZORFome2/el_merged.fasta", 'fasta')
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()

#part that determines SL in transcript
for sequence in intranscripts:
	name = sequence.name
	seq = sequence.seq
	if name in names:
		#print(name, len(seq))
		monstr_dic[name][11] = seq
		monstr_dic[name][4] = len(seq)
		if mySL in sequence.seq[:30] or revSL in sequence.seq[-30:]:
			monstr_dic[name][6] = "+"
		else:
			monstr_dic[name][6] = "-"


#print("Now printing control table, number of items per dictionary key...")
#for sequence in monstr_dic.keys():
#	print(sequence, len(monstr_dic[sequence]))


#########################
#### Printing Output ####
#########################

print("Now writing results...")
outtable = 'FINALTABLE-' + filename + '.txt'
#df = pd.DataFrame(monstr_dic.values(), columns=["annotation", "short", "contig", "no.reads", "contig len", "avg_cov", "SL", "TMHMM", "PrediSi", "PredSL", "TP_follow_SP", "nt seq", "aa seq"])
#df.to_csv(outtable, sep="\t")

#solution without pandas
outtable = open('FINALTABLE-' + filename, 'w')
for name in names:
	value = monstr_dic[name]
	outtable.write("{}\t{}\n".format(value[2],value[6]))
outtable.close()

print("Writing done. You are welcome.")
