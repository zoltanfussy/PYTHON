import argparse
import os

#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='geneious output to process', default="batch")
args = parser.parse_args()
infileparam = args.infile


codechange = {"Ala": "A", "Arg": "R", "Asn": "N",
"Asp": "D", "Cys": "C", "Glu": "E", 
"Gln": "Q", "Gly": "G", "His": "H",
"Ile": "I", "Leu": "L", "Lys": "K",
"Met": "M", "Phe": "F", "Pro": "P",
"Ser": "S", "Thr": "T", "Trp": "W",
"Tyr": "Y", "Val": "V", "Sup": "W", "Undet": "X"}

if infileparam == "batch":
	infiles = os.listdir('.')
	infiles = [f for f in infiles if f.endswith("all_modules.txt")]
	#example of these files in bico/genomeseq/spadesfilt/draft genome/mito_proteome/*
else: 
	infiles = [infileparam]

for infile in infiles:
	outfile = infile.split(".")[0] + ".gff3"
	c = 1
	with open(infile) as f, open(outfile, "w") as result:
		for line in f:
			data = line.split("\t")
			contig = data[0].replace(" ", "_").strip()
			if contig not in ("Sequence", "Name", "--------"):
				start = data[3].strip()
				stop = data[4].strip()
				gene = data[1].strip()
				#oneletter = codechange[tRNA]
				score = 0
				#no introns assumed
				if data[7].strip() == "forward":
					result.write("{0}\tgeneious\tgene\t{1}\t{2}\t{3}\t+\t.\tID={4}_{6};Name={4}_gene\n\
{0}\tgeneious\tCDS\t{1}\t{2}\t{3}\t+\t.\tID={4}{6};Name={4}\n\
{0}\tgeneious\texon\t{1}\t{2}\t{3}\t+\t.\tID={4}{6}_exon\n"
.format(contig, start, stop, score, gene, "", c))
					c += 1
				else:
					result.write("{0}\tgeneious\tgene\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6};Name={4}_gene\n\
{0}\tgeneious\tCDS\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6};Name={4}\n\
{0}\tgeneious\texon\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6}_exon\n"
.format(contig, start, stop, score, gene, "", c))
					c += 1

	print("Converting done.")


