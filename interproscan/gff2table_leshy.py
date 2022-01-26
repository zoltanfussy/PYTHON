import gzip,argparse,re

print("This script reformats interproscan gff3 files to annotation table")
print("usage: python gff2table_leshy.py <path/to/file>\n--------------------------------------------------------------------------")

#########################
####       Fx        ####
#########################

def read_panther_families(pantherfile="/mnt/mokosz/home/zoli/IPS/PANTHER16.0_HMM_classifications.gz"):
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
	line = string.split("\t")
	if line[2] == "polypeptide":
		writestring = ["", "", ""]
	elif line[1] in ["CDD", "Gene3D", "TIGRFAM", "PIRSF", "SFLD", "Coils", "ProDom", "SUPERFAMILY", "SMART"]:
		writestring = ["", "", ""]
	else:
		query, tool, annot = line[0], line[1], line[8]
		if ";" in annot:
			annots = annot.split(";")
			newannots = []
			if "signature_desc" in annot:
				has_name = True
			else:
				has_name = False
			for a in annots:
				if a.startswith("date") or a.startswith("Target") or a.startswith("ID"):
					continue
				elif a.startswith("signature_desc"):
					a = a.replace("signature_desc=", "")
					#put name in front
					newannots = [a] + newannots
				elif a.startswith("Name"):
					#if there is a valid signature description
					if has_name == True:
						#cannot have two names; this will be Pfam id
						continue
					#else if this is a panther family
					elif "PTHR" in a:
						name = panther_db.get(a.replace("Name=", ""), a)
						if name != a:
							newannots += [name,a.replace("Name=", "")]
						else:
							newannots += [a]
					else:
						#do not need other than name
						continue
			newannots = list(set(newannots))
		writestring = [query, tool, newannots]
	return(writestring)


def parse_fasta(seqname, sequence):
	if seqname.startswith(">match"):
		writestring = ""
	else:
		writestring = "{}\n{}\n".format(seqname, "\n".join(sequence))
	return(writestring)


def unique_list(inlist):
	unique_list = []
	for item in inlist:
		if item not in unique_list:
			unique_list.append(item)

	return unique_list


#########################
#### Read parameters ####
#########################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--interpro', help='InterProScan file to be processed', default='')
parser.add_argument('-b', '--blast', help='BLAST file to be processed', default='')
parser.add_argument('-f', '--fasta', help='FASTA file to be processed', default='')

args = parser.parse_args()

interprolist = args.interpro.split(",")
blastlist = args.blast.split(",")
fastalist = args.fasta.split(",")
uninformative = {"hypothetical", "putative", "N/A", "predicted:", "multispecies"}
contam_re = r"\{'(.*)'\}"

print("IPSparser: Files to be analyzed: {}".format(", ".join(interprolist)) )
print("BLASTparser: Files to be analyzed: {}".format(", ".join(interprolist)) )

global panther_db
panther_db = read_panther_families()

for file in interprolist:
	print("IPSparser: Opening file {}".format(file))
	filename = file.split(".")[0]
	with open(file) as f:
		print("IPSparser: Parsing annotations...")
		annotations = {}
		for line in f:
			line = line.strip()
			if line == "##FASTA":
				break

			#each annotation is one-line
			if line.startswith("##"):
				continue
			else:
				query, tool, newannots = parse_annotations(line)
				if tool.startswith("ProSite"):
					tool = "ProSite"
				if query == "":
					continue
				if query not in annotations:
					annotations[query] = {"ProSite": [], "Pfam": [], "PANTHER": [], "PRINTS": [], "BLAST": [], "BBH": "", "note": ""}
				try:
					annotations[query][tool].extend(newannots)
				except KeyError:
					continue
					#print("query", query)
					#print("tool", tool)

for file in blastlist:
	print("BLASTparser: Opening file {}".format(file))
	filename = file.split(".")[0]
	with open(file) as f:
		print("BLASTparser: Parsing annotations...")
		seen = set()
		for line in f:
			line = line.strip().split("\t")
			#each annotation is one-line
			if line[0].startswith("#"):
				continue
			else:
				query, newannots = line[0], line[-1] #modify here to exclude accession and organism
				newannots = newannots.split("[")[0].split()[1:]
				#if query == "":
				#	continue
				if query not in annotations:
					annotations[query] = {"ProSite": [], "Pfam": [], "PANTHER": [], "PRINTS": [], "BLAST": [], "BBH": line[1], "note": ""}
				elif query not in seen:
					seen.add(query)
					annotations[query].update({"BLAST": newannots, "BBH": line[1]})
				else:
					try:
						annotations[query]["BLAST"].extend(newannots)
					except KeyError:
						print("BLAST", "not in dictionary")
						continue

for file in fastalist:
	print("FASTAparser: Opening file {}".format(file))
	filename = file.split(".")[0]
	with open(file) as f:
		print("FASTAparser: Parsing annotations...")
		for line in f:
			#parse seq id
			if line.startswith(">"):
				query = line.split()[0].replace(">", "")
				newannots = line.replace(query, "")
				if "{" in newannots:
					newannots = re.search(contam_re, newannots).group(1)
				else:
					newannots = ""
				if query not in annotations:
					annotations[query] = {"ProSite": [], "Pfam": [], "PANTHER": [], "PRINTS": [], "BLAST": [], "BBH": "", "note": newannots}
				else:
					annotations[query].update({"note": newannots})


with open(filename + ".tsv", "w") as out:
	out.write("Protein\tProSite\tPfam\tPANTHER\tPRINTS\tBBH\tBLAST\tnote-contamination\n")
	for query in annotations:
		PSP = unique_list(annotations[query].get("ProSite", []))
		Pfam = unique_list(annotations[query].get("Pfam", []))
		PANTHER = unique_list(annotations[query].get("PANTHER", []))
		PRINTS = unique_list(annotations[query].get("PRINTS", []))
		BLAST = unique_list(annotations[query].get("BLAST", []))
		BLAST = [x for x in BLAST if x.lower() not in uninformative]
		BBH = annotations[query].get("BBH", "")
		NOTE = annotations[query].get("note", "")
		PSP = ";".join(PSP)
		Pfam = ";".join(Pfam)
		PANTHER = ";".join(PANTHER)
		PRINTS = ";".join(PRINTS)
		BLAST = " ".join(BLAST)
		out.write(f"{query}\t{PSP}\t{Pfam}\t{PANTHER}\t{PRINTS}\t{BBH}\t{BLAST}\t{NOTE}\n")

print("IPSparser: Finished.")		


