#this script renames accessions in final PASTA trees and alignments according to the name_translation file. A PASTA prefix needs to be given.
#optionally, unaligned input fasta can be renamed
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--prefix', help='PASTA prefix', required=True)
parser.add_argument('-f', '--fasta', help='FASTA file', default="default")

args = parser.parse_args()

prefix = args.prefix
if args.fasta != "default":
	fasta = args.fasta
else:
	fasta = prefix + ".fasta"

#first, name translation file is imported
print("PASTA prefix defined:%s" %(prefix))
intranslation = open(prefix + '_temp_name_translation.txt')
line = intranslation.read()
intranslation.close()

codes = line.split('\n')
nreads = len(codes) // 3

translations = {}
fullnames = {}
badchars = (",:.;()'")

#fill the reads list
for i in range(nreads):
	#this will repeat the loop i times, which equals the number of sequences in the dataset
	#NOTE TAKE only first 30 chars so the script is compatible with the prediction tree renamer:  codes[1][:30]
	fullname = codes[1]
	fullname = ''.join(c for c in fullname if c not in badchars)
	translations[codes[0]] = fullname
	fullnames[fullname] = ["sequence", "higher taxon", fullname]
	#removes first three lines
	codes = codes[3:]

infasta = SeqIO.parse(fasta, 'fasta')
notfound = 0
#notfoundlist = open(prefix + '_not-found.txt', 'w')

for seq in infasta:
	fullname = seq.name
	fullname = ''.join(c for c in fullname if c not in badchars)
	if fullname in fullnames:
		fullnames[fullname][0] = seq.seq
		try:
			highertaxon = seq.description.split(";")[2]
			taxon = seq.description.split(";")[-1]
			taxon = ''.join(c for c in taxon if c not in badchars)
			fullnames[fullname][1] = highertaxon
			fullnames[fullname][2] = taxon
		except:
			print("No higher taxon for accession: {}".format(seq.name))
			fullnames[fullname][1] = "unknown"
			fullnames[fullname][2] = fullname
	else:
		#notfoundlist.write(seq.name)
		notfound += 1
print("Full names NOT FOUND:" + str(notfound))
"""
for key in fullnames:
	print("{}\t{}".format(key, fullnames[key][2]))
"""

#reading the alignment to rename
print("Reading alignment to rename...")
inalign = open(prefix + '_temp_iteration_2_seq_alignment.txt')
line = inalign.read()
inalign.close()

lines = line.split('\n')
if len(lines) % 2 != 0:
	lines = lines[:-1]

outfile = open(prefix + '.alignment.fasta', 'w')
outfile2 = open(prefix + '.dataset.fasta', 'w')

for sequence in lines:
	if sequence.startswith('>'):
		safeheader = sequence[1:]
		header = translations[safeheader]
		newheader = (">{}_{}@{}\n".format(fullnames[header][2], header, fullnames[header][1]))
		outfile.write(newheader)
		outfile2.write(newheader)
	else:
		outfile.write(sequence + '\n')
		sequence = ''.join(c for c in sequence if c not in '-')
		outfile2.write(sequence + '\n')

outfile.close()
print("Alignment written to {}.alignment.fasta. Unaligned sequences are in {}.dataset.fasta. \n Reading tree to rename...".format(prefix, prefix))

strom = open(prefix + '_temp_iteration_2_tree.tre')
strom_line = strom.readline()

for key in translations:
	newkey = translations[key]
	newname = ("{}_{}@{}".format(fullnames[newkey][2], newkey, fullnames[newkey][1]))
	strom_line = strom_line.replace(key, newname)


with open(prefix + '.some.tre', 'w') as result:
    result.write(strom_line)

print("Tree written to {}.some.tree. Terminating.".format(prefix))