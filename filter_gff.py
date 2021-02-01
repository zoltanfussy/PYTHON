inputfile = "WW"

with open("{}.txt".format(inputfile)) as f:
	filterlist = set()
	for l in f:
		filterlist.add(l.strip())

with open("{}_Phaglo1.gff".format(inputfile), "w") as result, open("/Users/morpholino/OwnCloud/AndyLab/Phaeocystis/Phaeocystis globosa v2.3/interproscan/Phaglo1_1line_geneious.gff3") as f:
	file = f.read().split("##")
	result.write("##{}".format("##".join(file[1:4])))
	for item in file[4:]:
		if item.startswith("FASTA"):
			result.write("##FASTA\n")
			seqs = item.split(">")
			for seq in seqs:
				if any(x in seq for x in filterlist):
					result.write(">{}".format(seq))
		else:
			if any(x in item for x in filterlist):
				result.write("##{}".format(item))

with open("{}_Phaant1.gff".format(inputfile), "w") as result, open("/Users/morpholino/OwnCloud/AndyLab/Phaeocystis/Phaeocystis antarctica v2.2/interproscan/Phaant1_1line_geneious.gff3") as f:
	file = f.read().split("##")
	result.write("##{}".format("##".join(file[1:4])))
	for item in file[4:]:
		if item.startswith("FASTA"):
			result.write("##FASTA\n")
			seqs = item.split(">")
			for seq in seqs:
				if any(x in seq for x in filterlist):
					result.write(">{}".format(seq))
		else:
			if any(x in item for x in filterlist):
				result.write("##{}".format(item))
