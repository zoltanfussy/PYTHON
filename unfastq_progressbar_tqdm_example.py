from Bio import SeqIO
import os
from tqdm import tqdm
import subprocess
#handling both fasta and fastq inputs

homedir = "/Users/zoliq/"
#homedir = "/Volumes/zoliq data/"
wd = homedir + "Downloads/happyphatr"
os.chdir(wd)

files = [f for f in os.listdir('.') if f.endswith(".fasta")]

for f in tqdm(files, desc="Files"):
	CMD = "wc -l {}".format(f)
	output = subprocess.check_output(CMD, shell=True) #python read standard output of an issued command
	lines = int(str(output).split()[1])
	filetype = f.split(".")[-1]
	if filetype in ("fasta", "fas", "fst"):
		number_seqs = lines/2
	elif filetype == "fastq":
		number_seqs = lines/4
	inFile = SeqIO.parse(f, filetype)
	with open(f + '.whatever', 'w') as result:
		for sequence in tqdm(inFile, desc="sequences", total=number_seqs):
			name = sequence.description
			seq = sequence.seq
			result.write('>{}\n{}\n'.format(name,seq))

print("converting done")
