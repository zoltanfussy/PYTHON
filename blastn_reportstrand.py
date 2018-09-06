import os
from Bio import SeqIO

homedir = "/Users/zoliq/"
#homedir = "/Volumes/zoliq data/"
wd = homedir + "Downloads/happyphatr"
os.chdir(wd)
fasta = "contigs_ptwt.fasta"
fastadata = SeqIO.parse(fasta, "fasta")
blastout = "contigs_ptwt.out"
fwcontigs = "contigs_ptwt_allfw.fasta"
misscontigs = "contigs_ptwt_unmapped.fasta"
errorname = "contigs_ptwt_errors.txt"

##############################
###         FUNCS          ###
##############################

def decor(string):
	def wrap():
		print("===============")
		print(string)
		print("===============")
	return wrap


##############################
###          MAIN          ###
##############################

#first call blast:
print("Calling BLASTN")
os.system('blastn -db phatr3_ref.fasta -query {} -out {} -outfmt "6 qseqid qlen positive sstrand" \
	-max_hsps 1'.format(fasta, blastout))

print("Processing BLAST result")
best_hsps = {}
with open(blastout) as f:
	for l in f:
		data = l.split("\t")
		query = data[0]
		coverage = int(data[1]) / int(data[2])
		strand = data[3]
		if coverage > 0.5:
			best_hsps[query] = strand

missing = []
c = 0
with open(misscontigs, "w") as miss, \
open(fwcontigs, "w") as fw:
	for seq in fastadata:
		c += 1
		if c % 1000 == 0:
			print(c)
		if seq.name in best_hsps:
			if best_hsps[seq.name] == "plus":
				sequence = seq.seq
			else:
				sequence = seq.seq.reverse_complement()
			fw.write(">{}\n{}\n".format(seq.name, sequence))
		elif len(seq.seq) > 1000:
			missing.append(seq.name)
			miss.write(">{}\n{}\n".format(seq.name, seq.seq))

if len(missing) != 0:
	with open(errorname, "w") as errorfile:
		errorfile.write("Contigs missing a hit in db: {}\n".format(len(missing)))
		errorfile.write(", ".join(missing))
	decoration = decor("Errors occurred, please refer to: {}".format(errorname))
	decoration()