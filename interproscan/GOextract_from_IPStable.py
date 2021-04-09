import os

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
#wd = home + "Jankoviny/Tick_transcriptome/"
wd = home + "genomes/phaeocystis/Phaeocystis globosa v2.3/"

os.chdir(wd)

with open("GOtable.txt", "w") as result:
	data = open("IPStable.txt")
	for line in data:
		annot = line.split("\t")
		if len(annot) > 1:
			#print(annot)
			seq = annot[0]
			gos = annot[2]
			if gos != "":
				if " " in gos:
					gos = gos.replace(" ",",")
				result.write("{}\t{}\n".format(seq, gos))

c = 0
with open("GO-WEGO.txt", "w") as result:
	data = open("IPStable.txt")
	for line in data:
		c += 1
		if c % 1000 == 0:
			print(c)
		annot = line.split("\t")
		if len(annot) > 1:
			#print(annot)
			seq = annot[0]
			gos = annot[2]
			if gos != "":
				fingos = ""
				if " " in gos:
					gos = gos.split(" ")
					gos = ["GO:" + x for x in gos]
					gos.sort()
					finalgos = "\t".join(gos)
				else:
					finalgos = "GO:{}".format(gos)
				result.write("{}\t{}\n".format(seq, finalgos))

with open("GO-goatools.tsv", "w") as result:
	data = open("IPStable.txt")
	for line in data:
		c += 1
		if c % 1000 == 0:
			print(c)
		annot = line.split("\t")
		if len(annot) > 1:
			#print(annot)
			seq = annot[0]
			gos = annot[2]
			if gos != "":
				fingos = ""
				if " " in gos:
					gos = gos.split(" ")
					gos = ["GO:" + x for x in gos]
					gos.sort()
					finalgos = ";".join(gos)
				else:
					finalgos = "GO:{}".format(gos)
				result.write("{}    {}\n".format(seq, finalgos))