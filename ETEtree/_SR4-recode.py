import argparse
#or maybe 

# groups = ["AGNPST","CHWY","DEKQR","FILMV"] #SR4
# recode = {}
# for i,aas in enumerate(groups):
# 	for aa in aas:
# 		recode[aa] = str(i+1)

def find_format(suffix):
	accepted = ["fasta", "phylip", "clustal", "emboss", "nexus", "stockholm"]
	if suffix in accepted:
		return suffix
	elif suffix in ["fasta", "fna", "faa", "fas", "fa"]:
		return "fasta"
	elif suffix in ["ali", "aln"]:
		return "fasta" #probably
	elif suffix in ["phylip", "phy"]:
		return "phylip-relaxed"
	elif suffix in ["nexus", "nex"]:
		return "nexus"
	else:
		quit("unrecognized MSA format")


def _recode_alignment(filename, suffix):
	from Bio import AlignIO
	alignmentfile = AlignIO.read(filename, find_format(suffix))
	with open(filename.replace(suffix, "-SR4.aln"), "w") as result:
		for i,r in enumerate(alignmentfile):
			rname = r.id
			rseq = str(r.seq)
			for aa in recode:
				rseq = rseq.replace(aa, recode[aa])
			if i != len(alignmentfile) - 1:
				result.write(">{}\n{}\n".format(rname, rseq))
			else:
				result.write(">{}\n{}".format(rname, rseq))


def _recode_fasta(filename, suffix):
	from Bio import SeqIO
	sequencefile = [x for x in SeqIO.parse(filename, find_format(suffix))]
	with open(filename.replace("." + suffix, "-SR4.fasta"), "w") as result:
		for i,r in enumerate(sequencefile):
			rname = r.id
			rseq = str(r.seq)
			for aa in recode:
				rseq = rseq.replace(aa, recode[aa])
			if i != len(sequencefile) - 1:
				result.write(">{}\n{}\n".format(rname, rseq))
			else:
				result.write(">{}\n{}".format(rname, rseq))


def main():
#to parse arguments listed in the command after the script
	global recode
	recode = {'A': '1', 'G': '1', 'N': '1', 'P': '1', 'S': '1', 'T': '1', 
			'C': '2', 'H': '2', 'W': '2', 'Y': '2', 
			'D': '3', 'E': '3', 'K': '3', 'Q': '3', 'R': '3', 
			'F': '4', 'I': '4', 'L': '4', 'M': '4', 'V': '4',
			'-': '-', 'J': '-', 'X': '-'}
	parser = argparse.ArgumentParser(description='How to use argparse')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-f', '--fasta', help='Fasta input', default='')
	group.add_argument('-a', '--alignment', help='Alignment input', default='')

	args = parser.parse_args()
	if args.fasta != "":
		suffix = args.fasta.split(".")[-1]
		_recode_fasta(args.fasta, suffix)
	elif args.alignment != "":
		suffix = args.alignment.split(".")[-1]
		_recode_alignment(args.alignment, suffix)


if __name__ == '__main__':
	main()