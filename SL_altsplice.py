#skript na extrakci sekvenci s dino SL
print("This script filters an input fasta file based on a presence of a SL sequence anywhere within 50 nt from the beginning or end of a given sequence.")
print("Defining a --SL is optional, otherwise euglena SL is used.")
#print("Note that full length SL 'CCGTAGCCATTTTGGCTCAAG' seems too stringent, so 'TTTGGCTCAAG' is default here.")
#print("Suggested batch usage: for i in `ls *.fasta`;do python SLfilter.py -in $i -sl TTTGGCTCAAG;done")

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-in', '--infile', help='fasta to be processed', required=True)
parser.add_argument('-sl', '--SL', help='sequence of the spliced leader', default='TTTTTCG')

args = parser.parse_args()

inFile = SeqIO.parse(args.infile, 'fasta')
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()
out = open('output-' + args.infile, 'w')
counter = 0
bait2name_d = {}
baits = []


#THIS PART TO CREATE BAITS

noSLseqs = set()
for sequence in inFile:
	SLpresent = False
	if mySL in sequence.seq[:50]:
		SLpresent = True
		seq = sequence.seq
		#out.write('>{}\n{}\n'.format(name,bait))
	elif revSL in sequence.seq[-50:]:
		SLpresent = True
		seq = sequence.seq.reverse_complement()
		#out.write('>{}\n{}\n'.format(name,bait))
	else:
		noSLseqs.add(sequence.name)

	if SLpresent:
		counter += 1
		name = sequence.name
		bait  = seq.split(mySL)[1][:30]
		if len(bait) == 30:
			if bait not in bait2name_d:
				bait2name_d[bait] = [name]
				baits.append(bait)
			else:
				bait2name_d[bait].append(name)
with open("baits-el_merged2.txt", "w") as baitsfile:
	for bait in bait2name_d:
		baitsfile.write("{}@{}\n".format(bait, " ".join(bait2name_d[bait])))

		
with open("baits-el_merged2.txt") as f:
	baitsdata = f.read().split("\n")
	for line in baitsdata:
		if len(line) != 0:
			counter += 1
			bait = line.split("@")[0]
			name = line.split("@")[1]
			bait2name_d[bait] = [name]
			baits.append(bait)

print("Filtering done. Found {} sequences with SL.".format(counter))
print("Now to finding alternatively spliced 5` ends...")

inFile = SeqIO.parse(args.infile, 'fasta')
counter = 0
progress = 0
for sequence in inFile:
	progress += 1
	if progress % 1000 == 0:
		print(progress)
	# THIS WILL BE FASTER IF WE FILTER OUT CONTIGS THAT WE KNOW CONTAIN SL AT THEIR END
	try:
		if sequence.name in noSLseqs:
			for bait in baits:
				internal = mySL + bait
				if bait in sequence.seq and internal not in sequence.seq:
					out.write(">{} baited by {}@{}\n{}\n".format(sequence.name, bait, ' '.join(bait2name_d[bait]), sequence.seq))
					counter += 1
				elif bait in sequence.seq.reverse_complement() and internal not in sequence.seq.reverse_complement():
					out.write(">{} baited by {}@{}\n{}\n".format(sequence.name, bait, ' '.join(bait2name_d[bait]), sequence.seq.reverse_complement()))
					counter += 1
	except NameError:
		if 1 == 1:
			for bait in baits:
				internal = mySL + bait
				if bait in sequence.seq and internal not in sequence.seq:
					out.write(">{} baited by {}@{}\n{}\n".format(sequence.name, bait, ' '.join(bait2name_d[bait]), sequence.seq))
					counter += 1
				elif bait in sequence.seq.reverse_complement() and internal not in sequence.seq.reverse_complement():
					out.write(">{} baited by {}@{}\n{}\n".format(sequence.name, bait, ' '.join(bait2name_d[bait]), sequence.seq.reverse_complement()))
					counter += 1


print("{} possible alternatively spliced 5` ends found".format(counter))
out.close()

