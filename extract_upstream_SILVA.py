import argparse
import gzip
from Bio import SeqIO,Seq
from collections import Counter


def get_upstream(file, primer, reverse, upstream_length=2):
	upstream = []
	length = len(primer)
	primer = primer.replace("T", "U")
	matches = 0
	if reverse:
		primerseq = Seq.Seq(primer).reverse_complement()
		primer = str(primerseq)
		if file.endswith("gz"):
			with gzip.open(file, "rt") as handle:
				for seq in SeqIO.parse(handle, "fasta"):
					if primer in seq.seq:
						#rindex because I am looking at V9 at the end of the SSU
						matches += 1
						coordinate = str(seq.seq).rindex(primer) + length
						upstream.append(str(seq.seq[coordinate:coordinate + upstream_length].reverse_complement()).replace("U", "T"))
	else:
		if file.endswith("gz"):
			with gzip.open(file, "rt") as handle:
				for seq in SeqIO.parse(handle, "fasta"):
					if primer in seq.seq:
						#rindex because I am looking at V9 at the end of the SSU
						matches += 1
						coordinate = str(seq.seq).rindex(primer) - upstream_length
						upstream.append(str(seq.seq)[coordinate:coordinate + upstream_length].replace("U", "T"))
	c = Counter(upstream)
	#return c.most_common()
	freqs = []
	for i in "ATGC":
		for j in "ATGC":
			freqs.append((i+j, upstream.count(i+j)))
	print(matches, "matches")
	return freqs


def main():
	#to parse arguments listed in the command after the script
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-f', '--file', help='File to process', default="SILVA_138_SSURef_NR99_tax_silva.fasta.gz")
	parser.add_argument('-n', '--primername', help='Output prefix', default="")
	parser.add_argument('-p', '--primerseq', help='Primer sequence', default="")
	parser.add_argument('-r', '--reverse_complement', help='Reverse complement primer?', action='store_true') # to only look for presence/absence of this argument

	args = parser.parse_args()

	upstream = get_upstream(args.file, args.primerseq, args.reverse_complement)
	with open("_{}.tsv".format(args.primername), "wt") as result:
		for i, j in upstream:
			result.write("{}: {}\n".format(i, j))

if __name__ == "__main__":
	#main()
	pass


names = ["V9FW", "V9eug", "V9meta1", "V9meta2", "V9RV1", "V9RV2"]
primers = {"V9FW": "TTTGTACACACCGCCC", "V9eug": "ACCTTGTTACGACTTTTGC", "V9meta1": "GTTACGACTTCTCGTTCCT", "V9meta2": "GTTACGACTTCTCCTTCCT", "V9RV1": "CCTTCCGCAGGTTCACCTAC", "V9RV2": "CCTTCTGCAGGTTCACCTAC"}
reverse = [False, True, True, True, True, True]

for i, name in enumerate(names):
	print(name)
	upstream = get_upstream("SILVA_138_SSURef_NR99_tax_silva.fasta.gz", primers[name], reverse[i])
	with open("_{}.tsv".format(name), "wt") as result:
		for i, j in upstream:
			result.write("{}: {}\n".format(i, j))
