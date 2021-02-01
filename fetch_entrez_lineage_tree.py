import os
from Bio import SeqIO
import argparse
import re
import csv

from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()

#########################
#### Read parameters ####
#########################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be analyzed', default="batch")
parser.add_argument('-d', '--directory', help='Change working directory', default='.')

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
	defdir = "DolezalLab/SecY/"
	wd = home + defdir
	os.chdir(wd)
else:
	os.chdir(args.directory)


###############################
#### Open and parse inputs ####
###############################

allowed = ("fasta", "fas", "fst", "phy", "phylip")
if args.infile == "batch":
	infilelist = [x for x in os.listdir(".") if x.split(".")[-1] in allowed]
	infilelist = [x for x in infilelist if not x.startswith("safe")] #bc these have been created by a previous run
	infilelist = [x for x in infilelist if not x.startswith("trim")] #bc these have been created by a previous run
elif args.infile.split(".")[-1] in allowed:
	infilelist = [args.infile]
else:
	quit("file type not recognized - is it fasta/fas/fst or phy/phylip?")

print("FETCH_LINEAGES: Files to be analyzed: " + ", ".join(infilelist))

LINEAGEFILE = home + "progs/PYTHON/fetch_lineages.tsv"
with open(LINEAGEFILE, 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	tax2lineage = {r[0]: r[1] for r in reader}
missing = []

for file in infilelist:
	if file.split(".")[-1] in ("fasta", "fas", "fst"):
		data = SeqIO.parse(file, "fasta")
		for seq in data:
			genus = seq.name.split("_")[0]
			if genus not in tax2lineage:
				missing.append(genus)
			else:
				print("FETCH_LINEAGES included:", genus, tax2lineage[genus])
	elif file.split(".")[-1] in ("txt", "tsv"):
		with open(file) as f:
			for l in f:
				line = l.strip()
				genus = line.split()[0]
				if genus not in tax2lineage:
					 missing.append(genus)

name2taxid = ncbi.get_name_translator(missing)
#TESTING PURPOSES ONLY:
#missing = ['Homo', 'Aspergillus', 'Haloquadratum']
#taxid2name = ncbi.get_taxid_translator([9606, 9443])
#rankofnode = ncbi.get_rank([9606, 9443])

for genus in name2taxid:
	#retrieve at least one species:
	descendants = ncbi.get_descendant_taxa(genus)
	lineage = ncbi.get_lineage(descendants[0])[2:7]
	names = ncbi.get_taxid_translator(lineage)
	rank = [names[taxid] for taxid in lineage]
	if "Eukaryota" in rank:
		rank.remove("Eukaryota")
	print("{}\t{}".format(genus, "_".join(rank)))

