from Bio import SeqIO
import os

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "genomes/euglena longa/ZORFome2"
os.chdir(wd)

exclude = set()
with open("../submission/exclude.txt") as exclusionfile:
	for l in exclusionfile:
		exclude.add(l.split()[0])

infile = SeqIO.parse("el_merged.fasta", "fasta")
n,m,x=0,0,0
positives,negatives = set(),set()
with open("SLpos.txt", "w") as pos, open("SLneg.txt", "w") as neg: 
	for sequence in infile:
		if sequence.name not in exclude:
			if 'TTTTTCG' in sequence.seq[:35]:
				n += 1
				pos.write(sequence.name + "\n")
				positives.add(sequence.name)
			elif 'CGAAAAA' in sequence.seq[-35:]:
				m += 1
				pos.write(sequence.name + "\n")
				positives.add(sequence.name)
			else:
				x += 1
				neg.write(sequence.name + "\n")
				negatives.add(sequence.name)

print("forward: {}, reverse: {}, sum(SL): {}".format(n, m, n+m))
print("no SL: {}".format(x))

# all_annotated = set()
# #with open("GO-WEGO_clusters.txt") as infile, \
# with open("./interproscan/amino.tsv") as infile, \
# open("SLpos_GO.txt", "w") as pos, open("SLneg_GO.txt", "w") as neg:
# 	pos.write("!SeqID\tGO:GOterms\n")
# 	neg.write("!SeqID\tGO:GOterms\n")
# 	for l in infile:
# 		if l == "":
# 			continue
# 		name = l.split()[0]
# 		#the .tsv output is much more redundant, so we keep all the alternative orfs
# 		"""if name not in all_annotated:
# 			all_annotated.add(name)
# 			newname = name
# 		else:
# 			newname = name + "_1"
# 			if newname in all_annotated:
# 				print("warning, contig present more than twice: " + name)
# 			all_annotated.add(newname)
# 		if name in positives:
# 			pos.write(l.replace(name, newname))
# 		elif name in negatives:
# 			neg.write(l.replace(name, newname))		"""

# 		#this part for the .tsv output of IPS on 6-frame translated datasets
# 		transcriptname = name[:-2]
# 		all_annotated.add(transcriptname)
# 		if transcriptname in positives:
# 			pos.write(l)
# 		elif transcriptname in negatives:
# 			neg.write(l)

# 	for seq in positives - all_annotated:
# 		pos.write(seq + "\t\n")
# 	for seq in negatives - all_annotated:
# 		neg.write(seq + "\t\n")

print("analysis done")