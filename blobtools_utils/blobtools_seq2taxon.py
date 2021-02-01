from Bio import SeqIO
import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
defdir = "phaeocystis/outgroup/"
wd = home + defdir
os.chdir(wd)

"""
#for multiple taxon-specific files in a folder
files = os.listdir('.')
files = [f for f in files if f.endswith(".pep") and os.path.isfile(f)]

with open("excavata.accession2taxid", "w") as result:
	for file in files:
		c = 0
		taxname = [file.split(".")[0].replace("_"," ")]
		print(taxname)
		if taxname == ["Aduncisulcus paluster"]:
			taxid = 6666666
		#elif taxname == ["Aduncisulcus paluster"]:
		#	taxid = 666666
		else:
			taxid = ncbi.get_name_translator(taxname)[taxname[0]][0]
		print(taxid)
		for sequence in SeqIO.parse(file, "fasta"):
			c += 1
			result.write("{0}\t{0}\t{1}\t000000000\n".format(sequence.name, taxid))
			if c % 1000 == 0:
				print(c)
"""
#for a single compiled file
taxacodes = {"Emiliania-huxleyi": "Emiliania huxleyi", "Gephyrocapsa-oceanica": "Gephyrocapsa oceanica", 
			"Isochrysis-galbana": "Isochrysis galbana", "Isochrysis-sp": "Isochrysis sp. CCMP1244", 
			"Prymnesium-parvum": "Prymnesium parvum", "Pleurochrysis-carterae": "Pleurochrysis carterae CCMP645",
			"Prymnesium-polylepis": "Prymnesium polylepis", "Chrysochromulina": "Chrysochromulina tobin"}
taxids = {"Emiliania-huxleyi": 2903, "Gephyrocapsa-oceanica": 38817, "Isochrysis-galbana": 37099, "Isochrysis-sp": 40639, 
			"Prymnesium-parvum": 97485, "Pleurochrysis-carterae": 13221, "Prymnesium-polylepis": 72548, "Chrysochromulina": 1460289}

file = "hapto_clust.fasta" #fasta with species prefixes
with open("haptophyta.accession2taxid", "w") as result:
	for seq in SeqIO.parse(file, "fasta"):
		c = 0
		taxcode = seq.name.split("_")[0]
		#print(taxcode)
		taxname = taxacodes[taxcode]
		#print(taxname)
		if taxcode in taxids:
			taxid = taxids[taxcode]
		#elif taxname == ["Aduncisulcus paluster"]:
		#	taxid = 666666
		else:
			print(taxname)
			taxdic = ncbi.get_name_translator([taxname])
			taxid = taxdic[taxname][0]
			print(taxid)
		c += 1
		result.write("{0}\t{1}\t000000000\n".format(seq.name, taxid))
		if c % 1000 == 0:
			print(c)