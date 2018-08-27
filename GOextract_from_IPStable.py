import os

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "Jankoviny/Tick_transcriptome/"

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
					finalgos = "\t".join(gos)
				else:
					finalgos = "GO:{}".format(gos)
				result.write("{}\t{}\n".format(seq, finalgos))