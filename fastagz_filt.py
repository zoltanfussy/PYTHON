import argparse, gzip, time
from Bio import SeqIO

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print("starting {}".format(current_time))

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='FASTA to process', required=True)
parser.add_argument('-f', '--filterfile', help='Filter to apply', default="")
args = parser.parse_args()
infile = args.infile
if args.filterfile:
	filtset = {x.strip() for x in open(args.filterfile).readlines()}
	print(len(filtset), "seqs to keep")
else:
	filtset = set()
	print("No filter applied")

with gzip.open(infile, mode='rt') as handle, open(infile.replace(".gz", ""), "w") as subset:
	if filtset:
		print("filtering fasta...")
		for seq in SeqIO.parse(handle, "fasta"):
			if seq.name in filtset:
				subset.write(">{}\n{}\n".format(seq.name, seq.seq))
				filtset.remove(seq.name)
	else:
		print("writing fasta...")
		for seq in SeqIO.parse(handle, "fasta"):
			subset.write(">{}\n{}\n".format(seq.name, seq.seq))

if args.filterfile:
	print(len(filtset), "seqs not in original data")
