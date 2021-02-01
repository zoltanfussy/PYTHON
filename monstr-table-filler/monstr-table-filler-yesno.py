print("This script produces a table of transcript/protein model characteristics of putative Euglena longa plastid-targeted proteins.")
print("(SL anywhere within 30 nt from the beginning or end of a transcript, sequence length, targeting algorithm predictions)")
print("Defining a --SL is optional, otherwise Euglena gracilis/longa SL is used.")
print("Note that full length SL 'ACTTTCTGAGTGTCTATTTTTTTTCG' seems too stringent, so 'TTTTTCG' is default here.")
print("Suggested batch usage: for i in `ls *.fas`;do python monstr-table-filler.py -in $i -sl TTTTTCG;done")

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

#### Functions ####
###################

def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via raw_input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
		It must be "yes" (the default), "no" or None (meaning
		an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True,
			 "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '{}'" .format(default))

	while True:
		sys.stdout.write(question + prompt)
		if sys.platform in ["darwin", "win32"]:
			choice = input().lower()
		elif sys.platform.startswith("linux"):
			choice = raw_input().lower()
		else:
			print("unrecognized OS, check how to use raw input")
			choice = raw_input().lower()

		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")

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
defdir = "progs/PYTHON-DATA/monstr-table-filler/"
os.chdir(home + defdir)

inproteins = SeqIO.parse(args.infile, 'fasta')
names = []
outproteins = open('outprot-' + args.infile, 'w')
for sequence in inproteins:
	annotation = " ".join(sequence.description.split(" ")[1:])
	name = sequence.name
	name = name.split("_")[0]
	if "-reversed" in name:
		name = name.replace("-reversed", "")
	if name not in names:	
		names.append(name)
	else:
		name = sequence.name
	seq = sequence.seq
	monstr_dic[name] = 	[annotation, "abbrev", name, "3=reads", "4=length", "5=avg_cov", "6=SL", "7=TMHMM", "8=PrediSi", "9=PredSL", "no SP, TP N/D", "11=nt", seq]
	outproteins.write(">{}\n{}\n".format(name, seq))

#monstr_dic["title"] = ["annotation", "short", "contig", "no.reads", "contig len", "avg_cov", "SL", "TMHMM", "PrediSi", "PredSL", "TP_follow_SP", "nt seq", "aa seq"]
#monstr_dic["indices"] = [  "0",        "1",      "2",       "3",         "4",        "5",    "6",    "7",       "8",     "9",         "10",        "11",     "12"]

outproteins.close()

intranscripts = SeqIO.parse(home + "genomes/euglena longa/ZORFome2/el_merged.fasta", 'fasta')
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()

#part that determines SL in transcript
outtranscripts = open('outtrans-' + args.infile, 'w')
tobowtie = open('tobowtie.fasta', 'w')
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
		outtranscripts.write(">{}\n{}\n".format(name, seq))
		tobowtie.write(">{}\n{}\n".format(name, seq))

outtranscripts.close()
tobowtie.close()

#print("Now printing control table, number of items per dictionary key...")
#for sequence in monstr_dic.keys():
#	print(sequence, len(monstr_dic[sequence]))

query_yes_no("predisi.txt, predsl.txt and tmhmm.txt files ready?")

###########################################################################
#          read in predictions from PrediSi, PredSL, TMHMM:               #
###########################################################################
predisi = open("predisi.txt").readlines()
predsl = open("predsl.txt").readlines()
tmhmm = open("tmhmm.txt").readlines()

print("parsing PrediSi, PredSL, and TMHMM predictions and sequences into data dictionary")
for line in predisi:
	line = line.split("\t")
	#FASTA-ID	Score	Cleavage Position	Signal Peptide ?	Chart
	#[0]		[1]		[2]					[3]					[4]
	name = line[0].split(" ")[0] #name = line[0].split("_")[0].replace("-reversed", "")
	if name in names:
		if line[3] == "Y":
			monstr_dic[name][8] = line[5]
		else:
			monstr_dic[name][8] = "-"

for line in predsl:
	line = line.split("\t")
	#sequence id	cTP score	mTP score	SP score	prediction	cleavage site	thylakoid	lTP cl. site
	#[0]			[1]			[2]			[3]			[4]			[5]				[6]			[7]
	name = line[0].split(" ")[0] #name = line[0].split("_")[0].replace("-reversed", "")
	if name in names:
		name = line[0]
		if line[4] == "secreted":
			monstr_dic[name][9] = line[5]
		else:
			monstr_dic[name][9] = "-"


for line in tmhmm:
	line = line.split("\t")
	name = line[0]#.split("_")[0].replace("-reversed", "")
	#queryID	length	expAA 	First60	PredHel	Topology
	#[0]		[1]		[2]		[3]		[4]		[5]
	numhel = line[4].split('=')[1]
	monstr_dic[name][7] = numhel

###########################################################################
#                      transit peptide detection:                         #
###########################################################################

print("preparing output for MultiLoc2")

with open('cleaved-' + args.infile, 'a') as result:
	for protein in monstr_dic.keys():
	#	print(protein, monstr_dic[protein])
		try:
			SP_si = monstr_dic[protein][8]
			SP_sl = monstr_dic[protein][9]
			if SP_si == "-":
				if SP_sl != "-":
					SP_sl = int(SP_sl)
					SP_cleaved_sl = '>{}@PredSL\n{}\n'.format(seqname, monstr_dic[protein][12][SP_sl:])
					result.write(SP_cleaved_sl)
			elif SP_sl == "-":
				SP_si = int(SP_si)
				SP_cleaved_si = '>{}@PrediSi\n{}\n'.format(seqname, monstr_dic[protein][12][SP_si:])
				result.write(SP_cleaved_si)
			elif int(SP_si) == int(SP_sl):
				SP_si = int(SP_si)
				SP_sl = int(SP_sl)
				SP_cleaved = '>{}\n{}\n'.format(seqname, monstr_dic[protein][12][SP_si:])
				result.write(SP_cleaved)
			else:
				SP_si = int(SP_si)
				SP_sl = int(SP_sl)
				SP_cleaved_si = '>{}\n{}\n'.format(seqname, monstr_dic[protein][12][SP_si:])
				SP_cleaved_sl = '>{}@PredSL\n{}\n'.format(seqname, monstr_dic[protein][12][SP_sl:])
				result.write(SP_cleaved_si)
				result.write(SP_cleaved_sl)
		except ValueError as VE:
			with open('result_errors.txt', 'a') as errors:
				errors.write('SP cleavage error: {}\n{}\n\n'.format(monstr_dic[protein][5], str(VE)))
print(monstr_dic)
print("cleaved proteins written to {}... now running MultiLoc2 plant option to determine transit peptide probabilities.".format("cleaved-" + args.infile))

if query_yes_no("Do you wish to perform MultiLoc2 on cleaved sequences?") is True:
	os.system('python /home/manager/MultiLoc2/src/multiloc2_prediction.py -fasta={} -predictor=LowRes -origin=plant -result=pred-ML2-cleave.txt -output=simple'.format("cleaved-" + args.infile))
#use the following if your script failed, but you have made the new MultiLoc2 predictions
#if os.path.isfile("pred-ML2-cleave.txt") == False:
#	os.system('python /home/manager/MultiLoc2/src/multiloc2_prediction.py -fasta={} -predictor=LowRes -origin=plant -result=pred-ML2-cleave.txt -output=simple'.format("cleaved-" + args.infile))

print("MultiLoc2 finished, now adding predictions to data dictionary")

ML2 = open("pred-ML2-cleave.txt").read()
ML2C = ML2.split('\n')
alternativecleavage = open("alternativecleavage-" + filename + ".txt", "w")

for line in ML2C:
	if len(line.split('\t')) != 1:
		line = line.split('\t')
		preds = sorted(line[1:])
		name = line[0]
		#['chloroplast: ', 'cytoplasmic: ', 'mitochondrial: ', 'nuclear: ', 'secretory pathway: ']
		if "@PredSL" in name or "@PrediSi" in name:
			alternativecleavage.write("alternative cleavage for: {}\n".format("\t".join(line)))
			name = name.split("@")[0]
			monstr_dic[name][10] = "check alt-cleavage file"
		else:
			try:
				cTP = preds[0].split()[1]
				mTP = preds[2].split()[1]
				#print(cTP, mTP)
				monstr_dic[name][10] = "pt: {} (mt: {})".format(cTP, mTP)
			except ValueError as VE:
				with open('result_errors.txt', 'a') as errors:
					errors.write('ML2 processing error: {}\n{}\n\n'.format(line, str(VE)))
print("Transit peptide predictions made, check alternativecleavage" + filename + ".txt for ambiguos predictions.")
alternativecleavage.close()


if query_yes_no("Do you wish to perform bowtie?") is True:
#run bowtie using os:
	os.system('bowtie-build -f {} DATABASE'.format("outtrans-" + args.infile))
	os.system('bowtie -a -v 2 DATABASE -f E.longa_Raw_reads.fasta --suppress 1,2,4,5,6,7,8 > bowtie.txt')
#os.system('bowtie -a -v 2 DATABASE -f E.longa_Raw_reads.fasta --al mappedreads.txt --suppress 1,2,4,5,6,7,8 > bowtie.txt')

if os.path.isfile("bowtie.txt"):
	bowtie = open("bowtie.txt").read().split()
	bowtie_dic = {}
	for line in bowtie:
		line = line
		if line not in bowtie_dic.keys():
			bowtie_dic[line] = bowtie.count(line)
	for hit in bowtie_dic:
		monstr_dic[hit][3] = bowtie_dic[hit]
		monstr_dic[hit][5] = float(monstr_dic[hit][3]) * 80 / float(monstr_dic[hit][4])
else:
	print("Bowtie results not provided.")

#########################
#### Printing Output ####
#########################

print("Now writing results...")
outtable = 'FINALTABLE-' + filename + '.txt'
df = pd.DataFrame(monstr_dic.values(), columns=["annotation", "short", "contig", "no.reads", "contig len", "avg_cov", "SL", "TMHMM", "PrediSi", "PredSL", "TP_follow_SP", "nt seq", "aa seq"])
df.to_csv(outtable, sep="\t")
"""
#solution without pandas
outtable = ('FINALTABLE-' + filename, 'w')
for name, value in monstr_dic.items():
	outtable.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(value[0],value[1],value[2],value[3],value[4],value[5],value[6],value[7],value[8],value[9],value[10],value[11],value[12]))
outtable.close()
"""
print("Writing done. You are welcome.")