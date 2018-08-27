from Bio import SeqIO
import argparse
import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()


parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be analyzed', required=True)
args = parser.parse_args()

dataset = args.infile.split("_")[0]
infileslist = os.listdir('.')
infileslist = [f for f in infileslist if f.startswith(dataset)]

#load fasta
seq_d = {}
seq_set = set()
badchars = ("|@+,:;()'") #also []/
for infile in infileslist:
	indataset = SeqIO.parse(infile, 'fasta')
	with open("error.log", "a") as error, open(dataset + "_v2.fasta", "a") as outfasta:
		errors = False
		for sequence in indataset:
			fullname = sequence.description
			newname = []
			for c in fullname:
				if c in badchars:
					c = "_"
				newname.append(c)
			shortname = ''.join(newname)
			hightaxon = shortname.split("_")[-1]
			if hightaxon not in ["Acidobacteria", "Actinobacteria", "Alveolata", "Alphaproteobacteria", "Amoebozoa", "Apicomplexa", "Apusozoa", "Archaea", "Asgardgroup", "Bacteroidetes", "Chlamydiae", "Chlorarachniophyceae", "Chloroflexi", "Choanoflagellida", "Cryptophyta", "Cryptomonadales", "Cyanobacteria", "Deinococcus-Thermus", "Discosea", "Euryarchaeota", "Firmicutes", "Fungi", "Gemmatimonadetes", "Glaucocystophyta", "Gonyaulacales", "Gymnodiniales", "Haptophyceae", "Heterolobosea", "Ichthyosporea", "Mamiellophyceae", "Metazoa", "Mycetozoa", "Noctilucales", "Oxyrrhinales", "Peridiniales", "Planctomycetes", "Pyrenomonadales", "Rhizaria", "Rhodophyta", "Rhodothermaeota", "Spirochaetes", "Stramenopiles", "Synergistetes", "Thaumarchaeota", "Thermotogae", "Viridiplantae"]:
			# list of used taxons incomplete!
				try:
					descendants = ncbi.get_descendant_taxa(shortname.split("_")[0])
					lineage = ncbi.get_lineage(descendants[0])[2:]
					names = ncbi.get_taxid_translator(lineage)
					rank = [names[taxid] for taxid in lineage]
					if "Proteobacteria" in rank:
						hightaxon = rank[2]
						if hightaxon == "delta/epsilon subdivisions":
							hightaxon = rank[3]
						shortname = shortname.replace("_Proteobacteria", "")
					elif "FCB group" in rank:
						hightaxon = rank[3]
					elif "Opisthokonta" in rank:
						hightaxon = rank[2]
					elif "Terrabacteria group" in rank:
						hightaxon = rank[2]
					elif "Bacteria" in rank:
						hightaxon = rank[1]
					elif "Eukaryota" in rank:
						hightaxon = rank[1]
					else:
						hightaxon = "_".join(rank[0:2])
					print(shortname, hightaxon, rank[:4])
				except:
					hightaxon = ""
					errors = True
					error.write("Taxonomy could not be obtained for {}.\n".format(shortname))
				hightaxon = hightaxon.replace("/","-")
				shortname += "_{}".format(hightaxon)
			safeseq = str(sequence.seq).replace("*","")
			if shortname not in seq_d:
				if safeseq not in seq_set:
					seq_d[shortname] = (fullname, safeseq)
					outfasta.write(">{}\n{}\n".format(shortname, safeseq))
				else:
					errors = True
					error.write("duplicate sequence, skipping:\n{}\n{}\n".format(shortname, safeseq))
			else:
				errors = True
				error.write("seq name not unique, skipping:\n{0}\n{1}\n{0}\n{2}\n".format(shortname, safeseq, seq_d[shortname][1]))
			#print(">" + shortname + "\n" + seq_d[shortname])
if errors:
	print("Errors occurred during sequence read, please refer to error.log")

