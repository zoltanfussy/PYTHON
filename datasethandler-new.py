# Performs fasta header renaming, and calls mafft and trimal to align and trim datasets 
# automatically on all fasta files in the current folder
# Alternatively, set folder with -d

import os
from Bio import SeqIO,AlignIO
import argparse
import re

#########################
#### Read parameters ####
#########################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be analyzed', default="batch")
parser.add_argument('-a', '--aligner', help='Aligner', default='run_pasta.py')
parser.add_argument('-t', '--treemaker', help='Program for tree inference', default='none')
parser.add_argument('-d', '--directory', help='Change working directory', default='.')

args = parser.parse_args()

##################################
#### Create working directory ####
##################################

home = "/Users/zoliq/ownCloud/"
#home = "/Volumes/zoliq data/OwnCloud/"
os.chdir(home + "genomes/chromera/trees/cell_cycle/new_set")
os.chdir(args.directory)
for generation in range(1,15):
	if os.path.isdir("TREE" + str(generation)) == False:
		outdir = "TREE" + str(generation)
		break

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

print("PHYLOHANDLER: Files to be analyzed: " + ", ".join(infilelist))
print("PHYLOHANDLER: Data output to dir: " + outdir)

#remove bad datasets and bad chars from names
badchars = ("|@+,:;()'") #also []/
taxonpattern = r'\[(.+)\]'
errors = False
for file in infilelist:
	print("PHYLOHANDLER: Processing file: " + file)
	extension = file.split(".")[-1]
	filename = file.replace("." + extension, "")
	if extension in ("fasta", "fas", "fst"):
		indataset = SeqIO.parse(file, 'fasta')
	elif extension in ("phy", "phylip"):
		indataset = SeqIO.parse(file, 'phylip')
	else:
		continue
	#load fasta
	seq_d = {}
	seq_set = set()
	with open("error.log", "a") as error:
		for sequence in indataset:
			fullname = sequence.description
			newname = []
			for c in fullname:
				if c in badchars:
					c = "_"
				newname.append(c)
			fullrename = ''.join(newname) #this will be replaced by md5-based codes
			shortname = fullrename.split(" ")[0]
			taxonmatch = re.search(taxonpattern, fullrename)#.group(1)
			if taxonmatch:
				taxon = taxonmatch.group(1)
				if sequence.name.startswith(taxon):
					fullrename = fullrename.replace("[{}]".format(taxon), "")
					fullname = fullname.replace("[{}]".format(taxon), "")
			safeseq = str(sequence.seq).replace("*","")
			if shortname not in seq_d:
				if safeseq not in seq_set:
					seq_d[shortname] = (fullname, fullrename, safeseq)
				else:
					errors = True
					error.write("file:{}\tduplicate sequence, skipping:\n{}\n{}\n".format(file, shortname, safeseq))
			else:
				errors = True
				error.write("file:{0}\tseq ID not unique, skipping:\n{1}\n{2}\n{1}\n{3}\n".format(file, shortname, safeseq, seq_d[shortname][1]))
			#print(">" + shortname + "\n" + seq_d[shortname])
	if errors:
		print("PHYLOHANDLER: Errors occurred during sequence read, please refer to error.log")

	print("PHYLOHANDLER: Done loading sequences, now to writing safe_file...")
	with open("rename-{}.txt".format(filename), "w") as renaming, open("safe-{}.fasta".format(filename), "w") as safefile:
		for key,value in seq_d.items():
			renaming.write("{}\t{}\n".format(value[1][:50], value[0]))
			safefile.write(">{}\n{}\n".format(value[1][:50], value[2]))

	if args.aligner == "run_pasta.py":
		command = "{0} -d protein -i safe-{1}.fasta -j {1} -o {2}".format(args.aligner, filename, outdir)
	elif args.aligner == "mafft":
		command = "{0} --maxiterate 1000 --localpair --thread 4 {1} > safe-{1}.aln".format(args.aligner, filename)
	print("PHYLOHANDLER: issuing aligner\n" + command)
	os.system(command)

	#copy and rename PASTA alignment to current directory and issue trimal
	if args.aligner == "run_pasta.py":
		os.system("cp ./{1}/{0}.marker001.safe-{0}.aln ./safe-{0}.aln".format(filename, outdir))
		print("PHYLOHANDLER: issuing trimmer:\ntrimal -in safe-{0}.aln -out trim-{0}.aln -fasta -automated1".format(filename))
		os.system("trimal -in ./safe-{0}.aln -out trim-{0}.aln -fasta -gt 0.1".format(filename)) #-gappyout / -automated1 / -gt 0.3
	elif args.aligner == "mafft":
		print("PHYLOHANDLER: issuing trimmer:\ntrimal -in safe-{0}.aln -out trim-{0}.aln -fasta -automated1".format(filename))
		os.system("trimal -in ./safe-{0}.aln -out trim-{0}.aln -fasta -gt 0.1".format(filename)) #-gappyout / -automated1 / -gt 0.3		

	#open trimal-trimmed alignment for dumping any gaps-only sequences
	trimalignmentfile = AlignIO.read("trim-{0}.aln".format(filename), "fasta")
	outfile1, outfile2 = "trimfilt-{0}.fasta".format(filename), "trimfilt-{0}.phy".format(filename)
	#filter out any sequences that are gaps-only after trimming
	filtalignmentfile = [record for record in trimalignmentfile if record.seq.count("-") != len(record.seq) and len(record.seq) != 0]
	with open(outfile1, "w") as result:
		for index, r in enumerate(filtalignmentfile):
			#get rid of the trailing newline character at the end of file:
			if index != len(filtalignmentfile) - 1:
				result.write(">{}\n{}\n".format(r.description, r.seq))
			else:
				result.write(">{}\n{}".format(r.id, r.seq))
		#count number of remaining sequences and their length
		count, length = len(filtalignmentfile), len(r.seq)

	#convert trimfile to phylip format (phylobayes)
	if os.stat(outfile1).st_size > 0: #check for file size
		ffile = AlignIO.read(outfile1, "fasta")
		AlignIO.write(ffile, outfile2, "phylip-relaxed")
	else:
		print("#####\nWARNING: File {} has zero size\n#####".format(outfile1))

	print("PHYLOHANDLER: Automated trimming done. \n \
		Trimming produced a file with {} sequences of {} sites\n\n\
		####################################".format(count, length))

	if args.treemaker != "none":
		if args.treemaker in ["iqtree-omp", "iqtree"]:
			treecommand = "-m TEST -mset LG -nt AUTO -s trimfilt-{}.fasta".format(filename)
			print("PHYLOHANDLER: Issuing software for tree inference:\n{} {}".format(args.treemaker, treecommand))
		os.system("{} {}".format(args.treemaker, treecommand))

print("PHYLOHANDLER: All requested analyses finished. Hooray!")
if errors:
	print("PHYLOHANDLER: Errors occurred during sequence read, please refer to error.log")
