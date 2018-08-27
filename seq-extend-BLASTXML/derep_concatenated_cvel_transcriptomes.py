from Bio import SeqIO

fasta1 = SeqIO.parse("Chromera velia_MMETSP0290-trans.fasta", "fasta")
fasta2 = SeqIO.parse("Chromera velia_NCBI-Woehle.fasta", "fasta")
fastaout = open("Chromera velia_MMETSP0290-trans_NCBI-Woehle.fasta", "w")

uniqseq = {}
for seq in fasta2:
	uniqseq[seq.seq] = seq.name

for seq in fasta1:
	uniqseq[seq.seq] = seq.name

for item in uniqseq:
	fastaout.write(">{}\n{}\n".format(uniqseq[item], item))