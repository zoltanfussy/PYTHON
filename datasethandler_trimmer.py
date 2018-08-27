import argparse
import os
from Bio import SeqIO,AlignIO

print("Welcome to automated alignment trimmer.")
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be analyzed', required="True")

args = parser.parse_args()
os.chdir("/Users/zoliq/ownCloud/Jankoviny/Tick_transcriptome/TREE3")
suffix = args.infile.split(".")[-1]
if suffix in ("fasta", "fas", "fst"):
	filetype = "fasta"
elif suffix in ("phy", "phylip"):
	filetype = "phylip"
elif suffix in ("aln"):
	filetype = "fasta"
else:
	quit("file type not recognized - is it fasta/fas/fst or phy/phylip?")
filename = args.infile.split(".")[0]

#run trimal
os.system("trimal -in {0} -out trim-{1}.aln -{2} -automated1".format(args.infile, filename, filetype)) #-gappyout / -automated1 / -gt 0.3
print("issuing: trimal -in {0} -out trim-{1}.aln -{2} -automated1".format(args.infile, filename, filetype))

#open trimmed alignment for dumping any gaps-only sequences
trimalignmentfile = AlignIO.read("trim-{0}.aln".format(filename), "fasta")
outfile1, outfile2 = "trimfilt-{0}.fasta".format(filename), "trimfilt-{0}.phy".format(filename)
#filter out any sequences that are gaps-only after trimming
filtalignmentfile = [record for record in trimalignmentfile if record.seq.count("-") != len(record.seq) and len(record.seq) != 0]
with open(outfile1, "w") as result:
	for index, r in enumerate(filtalignmentfile):
		#get rid of the trailing newline character at the end of file:
		if index != len(filtalignmentfile) - 1:
			result.write(">{}\n{}\n".format(r.description, r.seq))
		else:
			result.write(">{}\n{}".format(r.id, r.seq))
	#count number of remaining sequences and their length
	count, length = len(filtalignmentfile), len(r.seq)

#convert trimfile to phylip format (phylobayes)
if os.stat(outfile1).st_size > 0: #check for file size
	ffile = AlignIO.read(outfile1, "fasta")
	AlignIO.write(ffile, outfile2, "phylip-relaxed")
else:
	print("#####\nWARNING: File {} has zero size\n#####".format(outfile1))

print("Automated trimming done. Hooray!\nTrimming produced a file with {} sequences of {} sites\n\n".format(count, length))