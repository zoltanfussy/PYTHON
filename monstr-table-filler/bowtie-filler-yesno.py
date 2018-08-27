print("This script produces a table of average coverage values.")
print("Suggested usage: python bowtie-filler-yesno.py -in infile -f fasta/table")
#python monstr-table-filler.py -in infile -f fasta/table for larger analysis

import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#os.system('source ~/.bash_profile')

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-in', '--infile', help='file to be processed', required=True)
parser.add_argument('-f', '--format', help='infile format', default='table')

args = parser.parse_args()

filename = args.infile.split(".")[0]

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
		raise ValueError("invalid default answer: '%s'" % default)

	while True:
		sys.stdout.write(question + prompt)
		choice = input().lower()
		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")

monstr_dic = {}
########################################################################
###      read in seq names to extract nt/aa models of interest:      ###
########################################################################
if args.format == "table":
	intable = open(args.infile).readlines()
	names = {}
	nameset = set()
	outtranscripts = open('tobowtie-' + filename + '.fasta', 'w')
	for index, sequence in enumerate(intable):
		try:
			name = sequence.split("\t")[0]
			annotation = sequence.split("\t")[1]
			abbrev = sequence.split("\t")[2]
			seq = sequence.split("\t")[3].strip()
		except IndexError:
			print(sequence)
			name = "IndexError"
		if name not in names:	
			names[index] = name
			if len(seq) != 0:
				outtranscripts.write(">{}\n{}\n".format(name, seq))
		else:
			names[index] = name
			print(name, " present more than once")
		monstr_dic[name] = 	[annotation, abbrev, name, "3=reads", len(seq), "5=avg_cov"]
		
elif args.format == "fasta":
	intranscripts = SeqIO.parse(args.infile, 'fasta')
	names = {}
	nameset = set()
	for index, sequence in enumerate(intranscripts):
		annotation = " ".join(sequence.description.split(" ")[1:])
		name = sequence.name
		name = name.split("_")[0]
		if "-reversed" in name:
			name = name.replace("-reversed", "")
		if name not in nameset:	
			nameset.add(name)
			names[index] = name
		else:
			print(name, " present more than once")
		seq = sequence.seq
		monstr_dic[name] = 	[annotation, "abbrev", name, "3=reads", len(seq), "5=avg_cov"]

else:
	quit("File format unrecognized.")
#monstr_dic["title"] = ["annotation", "short", "contig", "no.reads", "contig len", "avg_cov"]
#monstr_dic["indices"] = [  "0",		"1",	  "2",	   "3",		 "4",		"5"]

"""
print("Now printing control table, number of items per dictionary key...")
for sequence in monstr_dic.keys():
	print(sequence, len(monstr_dic[sequence]))

"""
###########################################################################
#								 bowtie:								 #
###########################################################################
if query_yes_no("Do you wish to perform bowtie?") is True:
#run bowtie using os:
	if args.format == "table":
		print("Issuing: bowtie-build -f {} DATABASE".format('tobowtie-' + filename + '.fasta'))
		os.system('bowtie-build -f {} DATABASE'.format('tobowtie-' + filename + '.fasta'))
	else:
		os.system('bowtie-build -f {} DATABASE'.format(filename + '.fasta'))
	print("Issuing: bowtie -a -v 2 DATABASE -f /path/to/elonga_allreads.fasta --al mappedreads.txt --suppress 1,2,4,5,6,7,8 > bowtie.txt")
	if os.path.isfile("/Volumes/zoliq data/Stiahnuté/elonga_allreads.fasta"):
		os.system('bowtie -a -v 2 DATABASE -f /Volumes/zoliq\ data/Stiahnuté/elonga_allreads.fasta --al mappedreads.txt --suppress 1,2,4,5,6,7,8 > bowtie.txt')
#		os.system('bowtie -a -v 2 DATABASE -f /Volumes/zoliq\ data/Stiahnuté/elonga_oldreads.fasta --al mappedreads-old.txt --suppress 1,2,4,5,6,7,8 > bowtie-old.txt')
	elif os.path.isfile("/Users/zoliq/Downloads/elonga_allreads.fasta"):
		os.system('bowtie -a -v 2 DATABASE -f ~/Downloads/elonga_allreads.fasta --al mappedreads.txt --suppress 1,2,4,5,6,7,8 > bowtie.txt')
#		os.system('bowtie -a -v 2 DATABASE -f ~/Downloads/elonga_oldreads.fasta --al mappedreads-old.txt --suppress 1,2,4,5,6,7,8 > bowtie-old.txt')
	else:
		print("======================\nWARNING: reads file not found!\n======================")
#for bowtie version 2 I used:
#os.system('bowtie2 -a -f --reorder -x DATABASE -U E.longa_Raw_reads.fasta --met-file ./bowtie.txt')
#note that max mismatches here is 1
#but this may not support --suppress reporting in the same sense as version 1

if os.path.isfile("bowtie.txt"):
	bowtie = open("bowtie.txt").read().split()
	bowtie_dic = {}
	for line in bowtie:
		#line = line
		if line not in bowtie_dic.keys():
			bowtie_dic[line] = bowtie.count(line)
	for hit in bowtie_dic:
		try:
			monstr_dic[hit][3] = bowtie_dic[hit]
			monstr_dic[hit][5] = float(monstr_dic[hit][3]) * 80 / float(monstr_dic[hit][4])
		except KeyError as KE:
			with open("bowtieerrorfile.txt", "a") as bowtieerrorfile:
				bowtieerrorfile.write("problems with {}\n".format(str(KE)))

else:
	print("======================\nWARNING: bowtie results not provided.\n======================")

print("Now writing results...")
"""
outtable = 'FINALTABLE-' + filename + '.txt'
df = pd.DataFrame(monstr_dic.values(), columns=["annotation", "short", "contig", "no.reads", "contig len", "avg_cov"])
df.to_csv(outtable, sep="\t")
"""
#solution without pandas
outtable = open('FINALTABLE-' + filename + '.tsv', 'w')
for index in names:
	name = names[index]
	outtable.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(monstr_dic[name][0],monstr_dic[name][1],monstr_dic[name][2],monstr_dic[name][3],monstr_dic[name][4],monstr_dic[name][5]))
outtable.close()

print("Writing done. Hooray.")
