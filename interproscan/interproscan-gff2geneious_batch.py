# Performs interproscan gff3 file reformatting for Geneious import

import os,gzip
import argparse

print("This script reformats interproscan gff3 files for ones suitable for Geneious import")
print("usage: python interproscan-gff2geneious_batch.py [-i infile.fasta/batch -d workdirectory]\n--------------------------------------------------------------------------")

#########################
####       Fx        ####
#########################

def read_panther_families(pantherfile="../PANTHER16.0_HMM_classifications.gz"):
	with gzip.open(pantherfile, "rt") as f:
		panther_db = {x.split("\t")[0]: x.split("\t")[1] for x in f}
	return panther_db


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
	line = string.strip().split("\t")
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
					#these are added InterProScan information
					pass
				elif a.startswith("signature_desc"):
					#signature description gives a good name, put that in front
					a = a.replace("signature_desc", "Name")
					newannots = [a] + newannots
				elif a.startswith("Name"):
					if is_newname == True:
						#if there is a valid Pfam signature description, we derive a name from there
						#since we cannot have two names in one line, this is the Pfam id
						a = a.replace("Name", "Id")
						newannots += [a]
					elif "PTHR" in a:
						#most likely this is a panther family
						#get name from panther_db based on panther code
						#then append Name and panther code if applicable
						a = a.replace("Name=", "")
						name = panther_db.get(a, a)
						if name != a:
							newannots += ["Name=" + name, "Id=" + a]
						else:
							newannots += ["Name=" + a]
					#whatever else, this is probably good as it is
					else:
						newannots += [a]
				else:
					newannots += [a]
			annot = ";".join(newannots)
			print(annot)
		else:
			print("; not in annotation?")

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
parser.add_argument('-d', '--directory', help='Change working directory', default='TEST')
parser.add_argument('-m', '--merge', help='Merge files', action='store_true')

args = parser.parse_args()

wd = "."
os.chdir(wd)
os.chdir(args.directory)

if args.infile == "batch":
	filelist = [x for x in os.listdir(".") if x.endswith("gff3")]
	filelist = [x for x in filelist if not x.endswith("geneious.gff3")]
else:
	filelist = [args.infile]

print("IPSparser: Files to be analyzed: {}".format(", ".join(filelist)) )

global panther_db
panther_db = read_panther_families()

if args.merge:
	seqbuffer = ""
	with open("merged_geneious.gff3", "w") as out:
		for i,file in enumerate(filelist):
			print("IPSparser: Opening file {}".format(file))
			with open(file) as f:
				f = f.readlines()
				print("IPSparser: Parsing annotations...")
				annotation = True #we start in the annotation section
				header = f[:3]
				if i == 0:
					out.write("".join(header))
				for line in f[3:]:
					line = line
					#now to information parsing
					if line.strip() == "##FASTA":
						#entering sequence section, turn off annotation-like parsing
						annotation = False			
					elif annotation == True:
						#each annotation is one-line
						if line.startswith("##"):
							writestring = line
						else:
							writestring = parse_annotations(line)
						out.write("{}".format(writestring))
					else:
						#in the sequence section, parse FASTA-like data
						if line.startswith(">match"):
							break
						else:
							seqbuffer += line
		#finishing with the last file, write fasta buffer
		print("IPSparser: Now for something completely different. Parsing sequences...")
		out.write("##FASTA\n")
		out.write(seqbuffer)
	
else:
	for file in filelist:
		print("IPSparser: Opening file {}".format(file))
		filename = file.split(".")[0]
		with open(file) as f, open(filename + "_geneious.gff3", "w") as out:
			print("IPSparser: Parsing annotations...")
			annotation = True #we start in the annotation section
			for line in f:
				line = line
				#now to information parsing
				if line.strip() == "##FASTA":
					#entering sequence section, turn off annotation-like parsing
					print("IPSparser: Now for something completely different. Parsing sequences...")
					out.write("##FASTA\n")
					annotation = False
				elif annotation == True:
					#each annotation is one-line
					if line.startswith("##"):
						writestring = line
					else:
						writestring = parse_annotations(line)
					out.write("{}".format(writestring))
				else:
					#in the sequence section, parse FASTA-like data
					if line.startswith(">match"):
						break
					else:
						out.write(line)

print("IPSparser: Finished.")		


