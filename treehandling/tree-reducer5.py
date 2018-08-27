import os
import re
import argparse
from Bio import SeqIO

#set working directory
os.chdir("/Users/zoliq/ownCloud/genomes/bico/genomeseq/spadesfilt/draft genome/trees")

print("This script removes unwanted branches from a dataset based on an input nexus tree.")
print("Mark the branches to be removed in the input tree by a colour (using FigTree). Please use basic colours or format found in your nex file. ")
print("usage: python tree-reducer.py -i eno.fasta -t testtree.nex -c all\n--------------------------------------------------------------------------")
#WIN:black = #-16777216, #000000; green = #-16737997, #009933

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-t', '--tree', help='Treefile for trimming', required=True)
parser.add_argument('-c', '--colour', help='Branch colour filter', default='all')

args = parser.parse_args()


if args.infile.split(".")[-1] in ("fasta", "fas", "fst"):
	indataset = SeqIO.parse(args.infile, 'fasta')
elif args.infile.split(".")[-1] in ("phy", "phylip"):
	indataset = SeqIO.parse(args.infile, 'phylip')
else:
	quit("file type not recognized - is it fasta/fas/fst or phy/phylip?")
intree = open(args.tree).read()
filtercolour = args.colour

basecolours = {'blue': '0000ff', 'brown': '996633', 'cyan': '00ffff', 'green': '00ff00', 'magenta': 'ff00ff', 'orange': 'ff8000', 'purple': '800080', 'red': 'ff0000' , 'yellow': 'ffff00'}
black = ['-16777216', '000000']
if filtercolour in basecolours:
	filtercolour = basecolours[filtercolour]
elif filtercolour == 'all':
	print("any colour accepted")
else:
	print("unknown filter, setting to 'user-defined'. taxa with unrecognized colour codes will be retained")

#load fasta
seq_d = {}
badchars = ("|@+,:;()'") #also []/
for sequence in indataset:
	shortname = sequence.description.replace(" ","_")
	newname = []
	for c in shortname:
		if c in badchars:
			c = "_"
		newname.append(c)
	shortname = ''.join(newname)
	#print(shortname)
	seq_d[shortname] = sequence.seq
	#print(">" + shortname + "\n" + seq_d[shortname])

print("done loading sequences")
#load taxa from tree
alltaxa = [] # a list of leaf labels with their color
pattern = r"\[&!color=.+"
skip = []
skippedc = 0
keptc = 0
treelines = intree.split('\n')[4:] #extract only leaf labels
for line in treelines:
	if line != ';':
		line = line.replace("'", "")
		line = line.replace("\t", "").replace("@","_").replace(" ","_").replace("|","_")
		#line = line.split("@")[0]
		#linecolour = line + colour
		#print(linecolour)
		alltaxa.append(line)
		if "[&!color" in line:
			colour = re.search(pattern, line).group()
			newcolour = line.split('[&!color=#')[1].replace("]", "")
			#print(taxon)
			if newcolour in black:
				print("black detected for %s, keeping this taxon" % (line))
				keptc += 1
			elif filtercolour == 'all':
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			elif newcolour == filtercolour:
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			else:
				print("unknown colour detected for %s, keeping this taxon" % (line))
				keptc += 1
			newtaxon = line.split('[&!color=#')[0]
			#print(newtaxon)
			alltaxa = [newtaxon if x==line else x for x in alltaxa] 
		else:
			#print(line)
			keptc += 1
			colour = ""

	else:
		break


print("done loading taxa, omitted taxa listed in omitted-%s" % (args.infile))
#write omitted taxa
with open('omitted-' + args.infile, 'w') as f:
	for taxon in skip:
		nohighertaxon = taxon.split("@")[0]
		nospacetaxon = taxon.replace(" ","_")
		if seq_d.get(taxon) != None:
			f.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
		elif seq_d.get(nohighertaxon) != None:
			f.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
		elif seq_d.get(nospacetaxon) != None:
			f.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
		else:
			print("!!!!!KEY ERROR for omitted", taxon)

print("writing filtered dataset...")
#write results
with open('filtered-' + args.infile, 'w') as out:
	for taxon in alltaxa:
		nohighertaxon = taxon.split("@")[0]
		nospacetaxon = taxon.replace(" ","_")
		if taxon not in skip:
			if seq_d.get(taxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
			elif seq_d.get(nohighertaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
			elif seq_d.get(nospacetaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
			else:
				print("!!!!!KEY ERROR for filtered", taxon)

print("WRITING DONE, \n{0} taxa kept in filtered-{2}, \n{1} taxa omitted in omitted-{2}".format(keptc, skippedc, args.infile))