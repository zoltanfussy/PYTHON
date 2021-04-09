print("Suggested usage: python integrate_SPpreds.py -i INFILE.fasta")

import argparse
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#os.system('source ~/.bash_profile')

#set working directory
homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud/"
wd = homedir + "Terka/Bakalarka/allseq_bestpredictors/predictions"
os.chdir(wd)
resultdir = homedir + "progs/plastNN-master"

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='fasta to be processed', default="../cvel_FINALSETinclEXT.fasta")

args = parser.parse_args()

filename = args.infile.split("/")[-1].split(".")[0]
#check for input prediction files
if os.path.isfile("predisi-{}.txt".format(filename)) == False:
	print("please provide PrediSi file named: predisi-{}.txt".format(filename))
	print("or issue 'java JSPP eukarya.smx {} predisi-{}.txt'".format(args.infile, filename))
	quit("FATAL ERROR: file missing")
if os.path.isfile("predsl-{}.txt".format(filename)) == False:
	print("please provide PrediSi file named: predsl-{}.txt".format(filename))
	print("or issue 'perl PredSL.pl {} 0 > predsl-{}.txt'".format(args.infile, filename))
	quit("FATAL ERROR: file missing")
if os.path.isfile("asafind-{}.txt".format(filename)) == False:
	print("please provide ASAfind file named: asafind-{}.txt".format(filename))
	print("or issue 'python2 ASAfind.py -f {} -p signalp-{}.txt -s'".format(args.infile, filename))
	quit("FATAL ERROR: file missing")

print("All required files found, proceeding...")
monstr_dic = {}
########################################################################
#		read in the list of positive/negative reference seqs		   #
########################################################################

refslist = open("refslist-{}.txt".format(filename)).readlines()
refpt = [x.split("\t")[0] for x in refslist if x.split("\t")[1].startswith("BIPARTITE")]
refother = [x.split("\t")[0] for x in refslist if not x.split("\t")[1].startswith("BIPARTITE")]
print("{} positive reference seqs".format(len(refpt)))

########################################################################
#		read in seq names to extract nt/aa models of interest:		   #
########################################################################
inproteins = SeqIO.parse(args.infile, 'fasta')
print("Reading sequence data...")
names = []
badchars = ("JOBUXZ")

for sequence in inproteins:
	annotation = " ".join(sequence.description.split(" ")[1:])
	name = sequence.name
	if "-reversed" in name:
		name = name.replace("-reversed", "")
	if name not in names:	
		names.append(name)
	seq = ''.join(c for c in str(sequence.seq) if c not in badchars)
	monstr_dic[name] = 	{"Annotation": annotation, "Name": name, "Seq": seq}

#print("Now printing control table, number of items per dictionary key...")
#for sequence in monstr_dic.keys():
#	print(sequence, len(monstr_dic[sequence]))

########################################################################
#		  read in predictions from PrediSi, PredSL, TMHMM:			   #
########################################################################
predisi = open("predisi-{}.txt".format(filename)).readlines()
predsl = open("predsl-{}.txt".format(filename)).readlines()
asafind = open("asafind-{}.txt".format(filename)).readlines()

report_errors = False
print("Parsing signal peptide predictions and sequences into data dictionary...")
newnames = []
for line in predisi:
	line = line.split("\t")
	#FASTA-ID	Cleavage Position	Signal Peptide ?	Score
	#[0]		[1]					[2]					[3]
	if line[0] in names:
		name = line[0]
		if line[2] == "Y" and int(line[1]) < 101: #sometimes we see signal very far in the protein
			monstr_dic[name].update({"PrediSi site": int(line[1]), "PrediSi SPscore": float(line[3])})
		else:
			monstr_dic[name].update({"PrediSi site": "-", "PrediSi SPscore": float(line[3])})
	else:
		newnames.add(name)
print("->PrediSi: done")

for line in predsl:
	line = line.split()
	#sequence id	mTP score	SP score	prediction	cleavage site
	#[0]			[1]			[2]			[3]			[4]			
	if line[0] in names:
		name = line[0]
		#line[3] for nonplant prediction - complex algae
		if line[3] == "secreted":
			monstr_dic[name].update({"PredSL site": int(line[4]), "PredSL SPscore": float(line[2])})
		else:
			monstr_dic[name].update({"PredSL site": "-", "PredSL SPscore": float(line[2])})
	else:
		newnames.add(name)
print("->PredSL: done")

for line in asafind:
	#Identifier	SignalP	ASAfind cleavage position	ASAfind/SignalP cleavage site offset	ASAfind 20aa transit score	ASAfind Prediction
	#[0]		[1]		[2]							[3]										[4]							[5]	
	
	if not line.startswith('#') and len(line) != 0:
		line = line.split('\t')
		name = line[0]
		if name in names:
			try:
				sp = int(line[2])
				asa = sp + int(line[3])
			except ValueError:
				sp = "-"
				asa = "-"
			monstr_dic[name].update({'SignalP site': sp, 'ASAfind site': asa})
	else:
		newnames.add(name)
print("->ASAfind/SignalP: done")

if len(newnames) > 0:
	print("ERROR, new sequence names appeared!")
	report_errors = True

######################################################################
#					  output data for PlastNN:						 #
######################################################################

os.chdir(resultdir)
print("preparing output for PlastNN")
resultdir = filename.split("_")[0]
if os.path.isdir(resultdir) == False:
	os.mkdir(resultdir)
	os.mkdir(resultdir + "/train_data")
	os.mkdir(resultdir + "/unlabeled_data")

# use created subfolder instead and name the files as required by plastNN
with open('{}/unlabeled_data/data.txt'.format(resultdir), 'w') as testfasta, \
	open('{}/unlabeled_data/tp.txt'.format(resultdir), 'w') as test, \
	open('{}/unlabeled_data/rna.txt'.format(resultdir), 'w') as testrna, \
	open('{}/train_data/negative.txt'.format(resultdir), 'w') as negfasta, \
	open('{}/train_data/neg_tp.txt'.format(resultdir), 'w') as neg_sp, \
	open('{}/train_data/positive.txt'.format(resultdir), 'w') as posfasta, \
	open('{}/train_data/pos_tp.txt'.format(resultdir), 'w') as pos_sp, \
	open('{}/train_data/rna.txt'.format(resultdir), 'w') as ref_rna, \
	open('{}/result_errors.txt'.format(resultdir), 'w') as errors:

	#RNA files contains RNA data from all, and separately for unlabeled, but do I need them?

	testfasta.write("id seq\n")
	test.write("id first_aa_of_tp\n") #only those with SP are listed!
	testrna.write("id hr5 hr10 hr15 hr20 hr25 hr30 hr35 hr40\n")
	negfasta.write("id seq\n")
	neg_sp.write("id first_aa_of_tp\n")
	posfasta.write("id seq\n")
	pos_sp.write("id first_aa_of_tp\n")
	ref_rna.write("id hr5 hr10 hr15 hr20 hr25 hr30 hr35 hr40\n")

	for protein in monstr_dic.keys():
		#print(protein, monstr_dic[protein])
		try:
			#print(protein)
			SP_si = monstr_dic[protein]["PrediSi site"]
			SP_sl = monstr_dic[protein]["PredSL site"]
			SP_sp = monstr_dic[protein].get("SignalP site", "nd")
			SP_asa = monstr_dic[protein].get("ASAfind site", "nd")
			if SP_sp == "nd": #if SignalP predictions missed this item
				errors.write("Sequence not analyzed by SignalP:\n>{}\n{}\n".format(protein, monstr_dic[protein]["Seq"]))
				SP_sp = "-"
				SP_asa = "-"
			cleavagesite = set([SP_si, SP_sl, SP_sp, SP_asa])
			#print(protein, cleavagesite)
			cleavagesite.discard("-")

			gene = protein.split("__")[0] #my datasets are named accordingly
			if len(cleavagesite) >= 1:
				if gene in refpt:
					#write to training pos file
					for x in cleavagesite: #to report all predicted cleavage sites
						pos_sp.write("{}_{} {}\n".format(protein, x, x))
						posfasta.write("{}_{} {}\n".format(protein, x, monstr_dic[protein]["Seq"]))
						ref_rna.write("{}_{} {}\n".format(protein, x, " ".join(8*["0"])))
				elif gene in refother:
					#write to training neg file
					for x in cleavagesite:
						neg_sp.write("{}_{} {}\n".format(protein, x, x))
						negfasta.write("{}_{} {}\n".format(protein, x, monstr_dic[protein]["Seq"]))
						ref_rna.write("{}_{} {}\n".format(protein, x, " ".join(8*["0"])))
				else:
					for x in cleavagesite:
						test.write("{}_{} {}\n".format(protein, x, x))
						testfasta.write("{}_{} {}\n".format(protein, x, monstr_dic[protein]["Seq"]))
						testrna.write("{}_{} {}\n".format(protein, x, " ".join(8*["0"])))
			#only those with SP are analyzed
			#else: 
			#	#write to result and outfasta file
			#	test.write("{} {}\n".format(protein, "-"))
			#	testfasta.write("{} {}\n".format(protein, monstr_dic[protein]["Seq"]))


				#print(",".join([str(x) for x in cleavagesite]))
		except ValueError as VE:
			errors.write('SP cleavage error: {}\n{}\n\n'.format(protein, str(VE)))
		except KeyError as KE:
			errors.write('SP cleavage error: {}\n{}\n\n'.format(protein, str(KE)))

	if report_errors:
		errors.write("New sequence names appeared!\n{}\n".format(", ".join(newnames)))

print("cleaved proteins written to plastNN-{}.fasta ".format(filename))

print("Writing done. Hooray.")
