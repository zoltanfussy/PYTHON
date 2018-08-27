import os
import re

#this script parses Interproscan gff3 to table to easily link sequences with their annotations
homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "Jankoviny/Tick_transcriptome/"

os.chdir(wd)

files = os.listdir('./interproscan/')
files = [f for f in files if f.endswith(".gff3")]
GOs_found = set()
GOs_annotation = {}
GOs_seqlist = {}
seq_table = {}
GOpattern = r'"GO:\d+'
IPSpattern = r'"InterPro:\w+'
FAMpattern = r'Name=\w+'
ANNOTpattern = r'signature_desc=.+'

c = 0
for f in files:
	c += 1
	print("Analyzing file {} of {}".format(c, len(files)))
	with open('interproscan/' + f) as file:
		data = file.read()
		annotations = data.split("##FASTA")[0]
		items = annotations.split("##sequence-region ")
		for i in items[1:]:
			query = i.split()[0]
			seqid = "_".join(query.split("_")[:-1])
			frame = query.split("_")[-1]
			GOs = set() #nejak pres RE? 
			for hit in re.findall(GOpattern, i):
				GO = hit.replace('"GO:', '')
				GOs.add(GO)
				GOs_found.add(GO)
				if GO in GOs_seqlist:
					GOs_seqlist[GO].add(seqid)
				else:
					GOs_seqlist[GO] = {seqid}
			IPSid = set()
			for hit in re.findall(IPSpattern, i):
				IPS = hit.replace('"InterPro:', '')
				IPSid.add(IPS)
			familyIDs = set()
			for hit in re.findall(FAMpattern, i):
				FAM = hit.replace('Name=', '')
				familyIDs.add(FAM)
			signature_desc = set()
			for hit in re.findall(ANNOTpattern, i):
				annot = hit.replace('signature_desc=', '').split(";")[0]
				signature_desc.add(annot)
			curdict = {"frame": frame, "GOs": GOs, "IPSid": IPSid, "family IDs": familyIDs, "signature_desc": signature_desc}
			if seqid not in seq_table:
				seq_table[seqid] = curdict
			else:
				#print("ERROR Duplicate seq entry: {} frames {} & {}".format(seqid, seq_table[seqid]["frame"], frame))
				newseqid = seqid + "_" + frame
				seq_table[newseqid] = curdict
			for line in i.split("\n"):
				if '"GO:' in line:
					match = re.search(GOpattern, line).group().replace('"GO:', '')
					signature_desc = re.search(ANNOTpattern, line).group().replace('signature_desc=', '').split(";")[0]
					if GO not in GOs_annotation:
						GOs_annotation[GO] = {signature_desc}
					else:
						GOs_annotation[GO].add(signature_desc)

with open("foundGOs.txt", "w") as result:
	print("Writing found GO terms to file...")
	result.write("GO\tfunction\tseq IDs\n")
	for k in sorted(GOs_seqlist.keys()):
		printids = " ".join(GOs_seqlist[k])
		function = GOs_annotation.get(k, "N/A")
		if isinstance(function, str):
			printfx = function
		elif isinstance(function, set):
			printfx = " ".join(function)
		else:
			print(function, type(function))
		result.write("{}\t{}\t{}\n".format(k, printfx, printids))

with open("IPStable.txt", "w") as result:
	print("Writing IPS table to file...")
	result.write("SeqID\tframe\tGOterms\tIPS ID\tFamily ID\tsignature description\n")
	for k in sorted(seq_table.keys()):
		printframe = seq_table[k]["frame"]
		printgo = seq_table[k]["GOs"]
		if len(printgo) == 0:
			printgo = ""
		else:
			printgo = " ".join(printgo)
		printips = seq_table[k]["IPSid"]
		if len(printips) == 0:
			printips = ""
		else:
			printips = " ".join(printips)
		printfamily = seq_table[k]["family IDs"]
		if len(printfamily) == 0:
			printfamily = ""
		else:
			printfamily = " ".join(printfamily)
		printsignature = seq_table[k]["signature_desc"]
		if len(printsignature) == 0:
			printsignature = ""
		else:
			printsignature = " ".join(printsignature)
		result.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(k, printframe, printgo, printips, printfamily, printsignature))

print("parsing finished")