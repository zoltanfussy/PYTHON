# Performs interproscan gff3 file reformatting for Geneious import

import os
import argparse

print("This script reformats interproscan gff3 files for ones suitable for Geneious import")
print("usage: python IPSparser.py [-i infile.fasta/batch -d workdirectory]\n--------------------------------------------------------------------------")

#########################
####       Fx        ####
#########################

def delbadchars(string):
	badchars = ("|+,:;()' ") #also []/@
	n = []
	for c in string:
		if c in badchars:
			c = "_"
		n.append(c)
	result = "".join(n)
	return result

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None

def parse_annotations(string):
	line = string.split("\t")
	if line[2] == "polypeptide":
		writestring = ""
	elif line[1] in ["CDD", "Gene3D", "TIGRFAM", "PIRSF", "SFLD", "Coils", "ProDom", "SUPERFAMILY", "SMART"]:
		writestring = ""
	else:
		target, program, database, left, right, evalue, direction, something, annot = line[0], "IPSparser", line[1], line[3], line[4], line[5], ".", line[7], line[8]
		if ";" in annot:
			annots = annot.split(";")
			newannots = []
			if "signature_desc" in annot:
				is_newname = True
			else:
				is_newname = False
			for a in annots:
				if a.startswith("date") or a.startswith("Target") or a.startswith("ID"):
					pass
				elif a.startswith("signature_desc"):
					a = a.replace("signature_desc", "Name")
					newannots = [a] + newannots
				elif a.startswith("Name"):
					if is_newname == True:
						a.replace("Name", "Id")
					else:
						newannots += [a]
				else:
					newannots += [a]
			annot = ";".join(newannots)

		writestring = "\t".join([target, program, database, left, right, evalue, direction, something, annot]) + "\n"
	return(writestring)

def parse_fasta(seqname, sequence):
	if seqname.startswith(">match"):
		writestring = ""
	else:
		writestring = "{}\n{}\n".format(seqname, "\n".join(sequence))
	return(writestring)

#########################
#### Read parameters ####
#########################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be analyzed', default="batch")
parser.add_argument('-d', '--directory', help='Change working directory', default='.')

args = parser.parse_args()

wd = "/mnt/mokosz/home/zoli/IPS/"
os.chdir(wd)
os.chdir(args.directory)

if args.infile == "batch":
	filelist = [x for x in os.listdir(".") if x.endswith("gff3")]
	filelist = [x for x in filelist if not x.endswith("geneious.gff3")]
else:
	filelist = [args.infile]

print("IPSparser: Files to be analyzed: {}".format(", ".join(filelist)) )

for file in filelist:
	print("IPSparser: Opening file {}".format(file))
	filename = file.split(".")[0]
	with open(file) as f, open(filename + "_geneious.gff3", "w") as out:
		print("IPSparser: Parsing annotations...")
		annotation = True #we start in the annotation section
		firstseq = True
		for line in f:
			line = line.strip()
			#now to information parsing
			if line == "##FASTA":
				annotation = False

			if annotation == True:
				#each annotation is one-line
				if line.startswith("##"):
					writestring = line + "\n"
				else:
					writestring = parse_annotations(line)
				out.write("{}".format(writestring))
			elif line == "##FASTA":
				#entering sequence section, turn off annotation-like parsing
				print("IPSparser: Now for something completely different. Parsing sequences...")
				out.write("##FASTA\n")
			else:
				#in the sequence section, parse FASTA-like data
				if line.startswith(">"):
					if not firstseq:
						#sequence loading done, new sequence header found
						writestring = parse_fasta(seqname, sequence)
						out.write(writestring)
					sequence = []
					seqname = line
				else:
					sequence.append(line)
					firstseq = False #now we have at least one

print("IPSparser: Finished.")		


