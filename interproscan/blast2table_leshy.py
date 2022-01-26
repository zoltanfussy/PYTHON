import sys,gzip

infile = sys.argv[1]

print("This script reformats blast output files to annotation table")
print("usage: python blast2table_leshy.py <path/to/file>\n--------------------------------------------------------------------------")

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

#########################
#### Read parameters ####
#########################

filelist = [infile]
uninformative = {"hypothetical", "putative"}

print("BLASTparser: Files to be analyzed: {}".format(", ".join(filelist)) )

for file in filelist:
	print("BLASTparser: Opening file {}".format(file))
	filename = file.split(".")[0]
	with open(file) as f, open(filename + ".blast.tsv", "w") as out:
		print("BLASTparser: Parsing annotations...")
		annotations = {}
		for line in f:
			line = line.strip().split("\t")
			#each annotation is one-line
			if line[0].startswith("#"):
				continue
			else:
				target, newannots = line[0], line[-1].split()
				newannots = [x for x in newannots if x.lower() not in uninformative]
				tool = "BLAST"
				#if target == "":
				#	continue
				if target not in annotations:
					annotations[target] = {"BLAST": []}
				try:
					annotations[target][tool].extend(newannots)
				except KeyError:
					continue
					#print("target", target)
					#print("tool", tool)
		out.write("Protein\tBLAST\n")
		for target in annotations:
			BLAST = list(set(annotations[target].get("BLAST", [])))
			BLAST = ";".join(BLAST)
			out.write(f"{target}\t{BLAST}\n")

print("BLASTparser: Finished.")		


