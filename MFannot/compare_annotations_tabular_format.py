import os

wd = ""

def overlaps(moduleA, moduleB):
	startA = moduleA[0]
	endA = moduleA[1]
	startB = moduleB[0]
	endB = moduleB[1]
	if startA <= startB <= endA:
		if startA <= endB <= endA:
			overlap = "embed"
		else:
			overlap = "overlap"
	elif startA <= endB <= endA:
		overlap = "overlap"
	elif startB <= startA <= endB and startB <= startB <= endB:
		overlap = "embed"
	else:
		overlap = "none"
	return overlap

def tblread(file):
	tbl_dic = {}
	with open(file) as f:
		for line in f:
			if len(line.split("\t")) > 1:
				line = line.split("\t")
				if line[0] != "":
					start = int(line[0].replace(">", "").replace("<", ""))
					stop = int(line[1].replace(">", "").replace("<", ""))
					if start < stop:
						direction = "forward"
					else:
						direction = "reverse"
					modulerange = (min(start, stop), max(start, stop))
				elif line[3] == "gene":
					modulename = line[4].strip()
					if modulename.startswith("trn"):
						modulename = modulename[:4]
					tbl_dic[modulename] = {"Range": [modulerange], "Database": ["MFannot"], "Strand": direction}
	return tbl_dic


def tsvread(file):
	tsv_dic = {}
	with open(file) as f:
		data = f.read().split("\n")
		taxon = data[0].split("\t")[0]
		for line in data:
			#print(line)
			if len(line.split("\t")) > 1:
				line = line.split("\t")
				direction = line[7]
				start = int(line[3].replace(">", "").replace("<", ""))
				stop = int(line[4].replace(">", "").replace("<", ""))
				modulerange = (min(start, stop), max(start, stop))
				if line[2] == "gene":
					modulename = line[1].replace(" gene", "").replace(" CDS", "")
					if modulename.startswith("trn"):
						modulename = modulename[:4]					
					tsv_dic[modulename] = {"Range": [modulerange], "Database": ["NCBI"], "Strand": direction}
	return tsv_dic, taxon

##############################
###          MAIN          ###
##############################

tsv_dic, taxon = tsvread("Vischeria_KU501221.tsv") # NCBI exported annotations
tbl_dic = tblread("Vischeria_KU501221.tbl") # MFannot file

for key,value in tbl_dic.items():
	if key in tsv_dic:
		if overlaps(tbl_dic[key]["Range"][0], tsv_dic[key]["Range"][0]) in ("embed", "overlap"):
			print("Module {} overlap".format(key))
			if not tsv_dic[key]["Strand"] == value["Strand"]:
				print("Something fishy, strand of {} {} vs {}".format(key, tsv_dic[key]["Strand"], value["Strand"]))
				tsv_dic[key]["Strand"] = "unclear"
		else:
			tsv_dic[key]["Range"].append(value["Range"][0])
			tsv_dic[key]["Database"].append(value["Database"][0])
			print("New module {} with range {}-{}".format(key, value["Range"][0][0], value["Range"][0][1]))
	else:
		tsv_dic[key] = value
		print("New module {} with range {}-{}".format(key, value["Range"][0][0], value["Range"][0][1]))

with open("Vischeria_KU501221_all_modules.txt", "w") as result:
	for key in sorted(tsv_dic.keys()):
		for r in tsv_dic[key]["Range"]:
			result.write("{}\t{}\tgene\t{}\t{}\t{}\t1\t{}\n".format(taxon, key, r[0], r[1], abs(r[1]-r[0]), tsv_dic[key]["Strand"]))
			#the "1\t" substring can be omitted, it comes from geneious annotation table and represents the number of intervals for a given annotation
