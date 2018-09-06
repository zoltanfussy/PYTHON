from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
import os

homedir = "/Users/zoliq/"
#homedir = "/Volumes/zoliq data/"
wd = homedir + "Downloads/happyphatr"
os.chdir(wd)

files = os.listdir('.')
files = [f for f in files if f.endswith(".taxified.out")]

for f in files:
	f_translated = f.replace("taxified", "taxid_translated")
	with open(f) as infile, open(f_translated, "w") as outfile:
		table = infile.read().split("\n")
		table = [x for x in table if len(x.split("\t")) > 1]
		print("sorting contigs based on BLAST hits...")
		c = 0
		for line in table:
			c += 1
			if c % 10000 == 0:
				print("{} queries processed".format(c))
			line = line.split("\t")
			if line[1] != "N/A": #this is never False
				taxid = line[1]
				if taxid != "no-hit" and taxid != "N/A":
					lineage = ncbi.get_lineage(taxid)[2:]
					names = ncbi.get_taxid_translator(lineage)
					rank = [names[taxid] for taxid in lineage]
					group = "_".join(rank[1:3])
					orgn = rank[-1]
					outfile.write("{}\t{}>{}\n".format(line[0], orgn, group))