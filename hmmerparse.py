import os,sys,argparse
from Bio import SeqIO
import gzip

#os.chdir("/Users/zoliq/ownCloud/progs/PYTHON/mydata")
"""
prefixes = [x for x in os.listdir("./alignments") if x.endswith(".ali")]
#prefixes = ["PF00134_rp75", "PF02984_rp75"]
for prefix in prefixes:
	#os.system("hmmbuild {0}.hmm {0}.ali".format(prefix))
	os.system("hmmsearch -o {0}.result.txt --notextw --tblout {0}.table.txt --cpu 2 {0}.hmm hmmdb.fa".format(prefix, prefix))
"""
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-d', '--database', help='input database', required=True)
parser.add_argument('-i', '--ignore_threshold', help='Ignore inclusion threshold in hmmer output', action="store_true")
args = parser.parse_args()

#read in seqs from DB fasta only if seqname in results:
files = [x for x in os.listdir(".") if x.endswith("table.txt")]
alltables = set()
for file in files:
	data = open(file).read().split()
	for item in data:
		alltables.add(item)

seq_d = {}
try:
	fastadb = args.database
except IndexError:
	fastadb = "EukProt_v2_renamed.faa"

if fastadb.endswith("gz"):
	with gzip.open(fastadb, "rt") as handle:
		for seq in SeqIO.parse(handle, "fasta"):
			if seq.name in alltables:
				seq_d[seq.name] = str(seq.seq).replace("*", "")			
else:
	fastafile = SeqIO.parse(fastadb, "fasta")
	for seq in fastafile:
		if seq.name in alltables:
			seq_d[seq.name] = str(seq.seq).replace("*", "")

#filter results and export fasta
for file in files:
	multidomain = set()
	written = set() #move multidomain, written into the for loop to have domain-specific files
	name = file.replace(".table.txt", "-dbhits.fasta")
	with open(file) as infile, open(name, "w") as result:
		for line in infile:
			if "Domain annotation for each sequence" in line:
				break
			if "inclusion threshold" in line:
				if args.ignore_threshold:
					continue
				else:
					break
			line = line.split()
			if len(line) < 9:
				continue
			if  line[0] != "#":
				try:
					evalue = float(line[0])
					seqname = line[8]
					#print(seqname)
					try:
						if "::" in seqname:
							#process in-house transcriptome seqIDs to be shorter
							writeseqname = "{}__{}_{}".format(seqname.split("::")[0], seqname.split("::")[-2], seqname.split("::")[-1])
						else:
							writeseqname = seqname
					except IndexError:
						writeseqname = seqname
						print("Could not parse seqname", writeseqname)
					if seqname in written:				
						pass
					elif evalue < 1: #which is always...
						result.write(">{}\n{}\n".format(writeseqname, seq_d[seqname]))
						written.add(seqname)
						multidomain.add(seqname)
					elif seqname in multidomain: #several domains found in the same seq
						result.write(">{}\n{}\n".format(writeseqname, seq_d[seqname]))
						written.add(seqname)
					else:
						multidomain.add(seqname)
				except ValueError:
					pass
					#this is not a score table line

