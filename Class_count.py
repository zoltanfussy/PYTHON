from Bio import SeqIO
import os

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "genomes/euglena longa/ZORFome2"
os.chdir(wd)

exclude = set()
with open("../submission/exclude.txt") as exclusionfile:
	for l in exclusionfile:
		exclude.add(l.split()[0])

classI,classII = set(),set()
for seq in SeqIO.parse("class_I.fa", "fasta"):
	if seq.name not in exclude:
		classI.add(seq.name)
for seq in SeqIO.parse("class_II.fa", "fasta"):
	if seq.name not in exclude:
		classII.add(seq.name)

with open("GO-WEGO_clusters.txt") as infile, \
open("ClassI_GO.txt", "w") as pos, open("ClassII_GO.txt", "w") as neg:
	pos.write("!SeqID\tGO:GOterms\n")
	neg.write("!SeqID\tGO:GOterms\n")
	for l in infile:
		if l == "":
			continue
		if l.split()[0] in classI:
			pos.write(l)
		elif l.split()[0] in classII:
			neg.write(l)