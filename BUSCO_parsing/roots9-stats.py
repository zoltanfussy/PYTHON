#!/usr/bin/env python
import os
from Bio import SeqIO

#this script actually needs a lot of files with BUSCO full results

####################################
###########  SET INPUTS  ###########
####################################
#read list of all buscos
buscolistfile = open("BUSCOlist.txt")
#optionally, provide a list of BUSCOs that will be filtered to a special table
FILTLIST = ["EOG0937019B", "EOG093701S0", "EOG0937036L", "EOG0937050C", "EOG093705YE", "EOG093705Z6", "EOG09370ARO", 
"EOG09370AV1", "EOG09370B7D", "EOG09370D13", "EOG09370F47", "EOG09370GC2", "EOG09370L23", "EOG09370LFI", "EOG09370M9E", 
"EOG09370MMK", "EOG09370OGU", "EOG09370OLY", "EOG09370P2W", "EOG09370PYJ", "EOG09370QJ1", "EOG09370RCQ", "EOG09370SYK", 
"EOG09370TCX", "EOG09370UEC", "EOG09370UO4", "EOG09370UR3", "EOG09370USJ", "EOG09370X1T", "EOG09370XFJ", "EOG09370YCN", 
"EOG09371088", "EOG09371BIR", "EOG093705DM", "EOG09370762", "EOG093707RF", "EOG09370BJE", "EOG09370CRF", "EOG09370FLD", 
"EOG09370IF5", "EOG09370U0A", "EOG09370W6N"] 
#to obtain the FILTLIST, I set filter on <1.2 copies per annotation, <0.1 missing, <0.1 frgm, and <0.16 dupe data in Excel
#then it is probably a good idea to remove short alignments - see root9-stats.xls
#duplicates have been overestimated, as datasets contain a low number of duplicate sequences that should have been removed prior to BUSCO analysis
FILTSPECFILE = open("FILTSPEC.txt")
#leave FILTSPEC.txt empty when you don`t want any species filtered


#PROCESS INPUTS:
BUSCO_LIST = buscolistfile.read().split("\n")
print("BUSCOs searched: ", len(BUSCO_LIST))
buscolistfile.close()

FILTSPEC = FILTSPECFILE.read().split("\n")
FILTSPECFILE.close()

#determine list of BUSCO report files in current folder
filenames = os.listdir('.')
filenames = [f for f in filenames if f.endswith(".fasta.tsv")]
#skip = ["unwanted_species"]
filecount = len(filenames)


####################################
#############  MAIN:  ##############
####################################
#create a dictionary of BUSCO stats across all datasets
BUSCO_DICT = {}
for item in BUSCO_LIST:
	BUSCO_DICT[item] = {"Complete": set(), "Duplicated": set(), "Fragmented": set(), "Missing": set(), "Copies": 0}

#create a list of all BUSCOs found across all datasets to enable extraction of FASTAs
SPECIES_DICT = {}
orthologdict = {}
ortholognames = set()

outfastalist = open("roots9buscos-list.txt", "w")
#go through BUSCO reports (stored in filenames) and create stats for BUSCOs/datasets
for file in filenames:
	dataset = file.split(".fasta")[0].replace("full_table_","")
	reportfile = open(file)
	report = reportfile.read().split("\n")
	reportfile.close()
	curdict = {}
	curset = set()
	for line in report:
		if not line.startswith("#") and len(line) > 0:
			line = line.split()
			BUSCOID = line[0]
			if len(line) != 2:
				fastaid = line[2]
				outfastalist.write(fastaid + "\n")
				ortholognames.add(fastaid)
				orthologdict[fastaid] = BUSCOID
			status = line[1]
			BUSCO_DICT[BUSCOID][status].add(dataset)
			BUSCO_DICT[BUSCOID]["Copies"] += 1
			if status in ["Complete"]:
				curdict[BUSCOID] = "@"
				curset.add(BUSCOID)
			elif status in ["Fragmented"]:
				curdict[BUSCOID] = "f"
				curset.add(BUSCOID)
			elif status in ["Duplicated"]:
				curdict[BUSCOID] = "2"
				curset.add(BUSCOID)
			else:
				curdict[BUSCOID] = "-"
	curdict["TOTAL"] = len(curset)
	SPECIES_DICT[dataset] = curdict
outfastalist.close()
print("BUSCOs in dataset: ", len(ortholognames))

outspeciestable = open("roots9-spectable.tsv", "w")
outspeciestable.write("species\ttotal\t{}\n".format("\t".join(BUSCO_LIST)))

outfilteredtable = open("roots9-filttable.tsv", "w")
outfilteredtable.write("species\t{}\tfull single-copy\ttotal (of {} filtered)\n".format("\t".join(FILTLIST), len(FILTLIST)))
#this part writes the table of presence/absence of BUSCOs per dataset 
for species in SPECIES_DICT:
	outspeciestable.write("\n{}\t{}\t".format(species, SPECIES_DICT[species]["TOTAL"]))
	#write to speciestable presence/absence of busco by busco in a tab-separated format
	for BUSCO in BUSCO_LIST:
		outspeciestable.write(SPECIES_DICT[species][BUSCO] + "\t")
	if FILTLIST != []:
		outfilteredtable.write("{}\t".format(species))
		filtcount = 0
		singlecopy = 0
		for BUSCO in FILTLIST:
			present = SPECIES_DICT[species][BUSCO]
			if present == "@":
				filtcount += 1
				singlecopy += 1
			elif present == "2":
				filtcount += 1
			elif present == "f":
				filtcount += 1
			outfilteredtable.write(present + "\t")
		outfilteredtable.write("{}\t{}\n".format(singlecopy, filtcount))
	else:
		outfilteredtable.write("ERROR: No BUSCO filter defined...")

"""
#this part prepares stats of BUSCOs as for their duplicity and fragmentation
outtable = open("roots9-stats.tsv", "w")
outtable.write("ID\tcomp\tfrgm\tnodupe\tdupe\tmiss\tcp pa\tcontrol\n")
for key in BUSCO_DICT:
	Complete = len(BUSCO_DICT[key]["Complete"])
	completefract = "%.3f" % (Complete/filecount)
	Fragmented = len(BUSCO_DICT[key]["Fragmented"])
	frgmfract = "%.3f" % (Fragmented/filecount)
	Nodupe = Complete + Fragmented
	nodupefract = "%.3f" % (Nodupe/filecount)
	Duplicated = len(BUSCO_DICT[key]["Duplicated"])
	dupefract = "%.3f" % (Duplicated/filecount)
	Missing = len(BUSCO_DICT[key]["Missing"])
	missfract = "%.3f" % (Missing/filecount)
	Copies = BUSCO_DICT[key]["Copies"]
	copiesfract = "%.2f" % (Copies/filecount)
	controlcount = Complete + Fragmented + Duplicated + Missing
	outtable.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, completefract, frgmfract, nodupefract, dupefract, missfract, copiesfract, controlcount))

outtable.close()
#to obtain the FILTLIST, I set filter on <1.2 copies per annotation, <0.1 missing, <0.1 frgm, and <0.125 dupe data in Excel
"""

#this part takes all datasets and extracts BUSCO sequences into a single file
#and then distributes sequences to individual files grouped by their BUSCO assignment after pooling
#filters may be applied

#modify path below according to the folder that contains the files
path = os.fspath('/Users/zoliq/ownCloud/genomes/sources/euks/')
#path = os.fspath('/Volumes/zoliq data/OwnCloud/genomes/sources/euks/')
datasets = os.listdir(path=path)
datasets = [f for f in datasets if f.endswith(".fasta")]
#print(datasets)
"""
#uncomment this section to export all sequences from all BUSCO datasets
#i use seqdict to dereplicate sequences
outfasta = open("roots9-buscosonly.fas", "w")
seqdict = {}
folder1 = "allsets-allbuscos/"
for file in datasets:
	fileaddress = path + file
	print("searching ", fileaddress)
	for sequence in SeqIO.parse(fileaddress, "fasta"):
		if sequence.name in ortholognames:
			ortholog = orthologdict[sequence.name]
			if sequence.seq not in seqdict:
				seqdict[sequence.seq] = ("{}_{}".format(ortholog, sequence.name))
				outfasta.write(">{}_{}\n{}\n".format(ortholog, sequence.name, sequence.seq))
				with open(folder + ortholog + ".fasta", "a") as orthologset:
					orthologset.write(">{}_{}\n{}\n".format(ortholog, sequence.name, sequence.seq))
			else:
				print("{} from group {} duplicated, same as {}. Check for transcriptome contamination?".format(sequence.name, ortholog, seqdict[sequence.seq]))
"""

#this section prepares filtered BUSCO fastas or dies if no filter is set by FILTSPEC.txt
if FILTSPEC == ['']:
	print("analysis finished")
	quit("quitting because no datasets to filter")
else:
	print("continuing with BUSCO/dataset filtering and sorting")
	#otherwise continue to filter BUSCOs

#anotherseqdict to remove duplicate sequences
anotherseqdict = {}
folder = "FILTSPEC-FILTBUSCOS/"
errorfile = open("root9-dupe-error.log", "w")
errorfile2 = open("root9-warning.log", "w")
counter = 0
for file in datasets:
	counter += 1
	if counter % 10 == 0:
		print(str(counter) + " files sorted")	
	transcriptome = file.split(".fasta")[0]
	if transcriptome in FILTSPEC:
		fileaddress = path + file

		for sequence in SeqIO.parse(fileaddress, "fasta"):
			if sequence.name in ortholognames and orthologdict[sequence.name] in FILTLIST:
				ortholog = orthologdict[sequence.name]
				if sequence.seq not in anotherseqdict:
					anotherseqdict[sequence.seq] = ("{}_{}".format(ortholog, sequence.name))

##########UNCOMMENT HERE TO WRITE THE FASTAS##################
				#	with open(folder + ortholog + ".fasta", "a") as orthologset:
				#		orthologset.write(">{}\n{}\n".format(sequence.name, str(sequence.seq).replace("*", "")))
				else:
					errorfile.write("{} = {}\n".format(sequence.name, anotherseqdict[sequence.seq]))
					first = sequence.name.split("_")[0]
					second = anotherseqdict[sequence.seq].split("_")[1]
					if first != second:
						print("organism does not match!")
						errorfile2.write("{} same as {} but organism does not match. Check for transcriptome contamination?\n".format(sequence.name, anotherseqdict[sequence.seq]))
					if not anotherseqdict[sequence.seq].startswith(ortholog):
						print("BUSCO does not match!")
						errorfile2.write("{} (BUSCO {}) same as {} but BUSCO does not match. Check for transcriptome contamination or fusion proteins?\n".format(sequence.name, ortholog, anotherseqdict[sequence.seq]))

print("see *-dupe-error.log and *-warning.log for found duplications")
print("sorting done")

