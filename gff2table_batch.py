import os
import re

#this script parses Interproscan gff3 to table to easily link sequences with their annotations
if os.path.isdir("/Users/morpholino/OwnCloud/"):
	homedir = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	homedir = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
defdir = "Jankoviny/Tick_transcriptome/"
defdir = "AndyLab/Phaeocystis/annotations_comparison"
wd = homedir + defdir

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
	print("Processing file {} of {}".format(c, len(files)))
	with open('interproscan/' + f) as file:
		data = file.read()
		annotations = data.split("##FASTA")[0]
		items = annotations.split("##sequence-region ")
		for i in items[1:]:
			query = i.split()[0]
			seqid = "_".join(query.split("_")[:-2]) #:-1 if only the last part >> frame is to be removed
			frame = query.split("_")[-1]
			#first find regex in the whole item:
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
			#find the range covered by annotations:
			lower, upper = [], [] 
			for line in i.split("\n"):
				columns = line.split("\t")
				#print(columns)
				if len(columns) == 1:
					pass
				elif columns[2] != "polypeptide": #this covers the whole sequence
					lower.append(int(columns[3]))
					upper.append(int(columns[4]))
			seqrange = "{}-{}".format(min(lower), max(upper))
			curdict = {"frame": frame, "range": seqrange, "GOs": GOs, "IPSid": IPSid, "family IDs": familyIDs, "signature_desc": signature_desc}
			if seqid not in seq_table:
				seq_table[seqid] = curdict
			else:
				#print("ERROR Duplicate seq entry: {} frames {} & {}".format(seqid, seq_table[seqid]["frame"], frame))
				#uncomment this part to create a new entry:
				#newseqid = seqid + "_" + frame
				#seq_table[newseqid] = curdict
				#uncomment this part to fuse:
				seq_table[seqid]["GOs"] |= GOs
				seq_table[seqid]["IPSid"] |= IPSid
				seq_table[seqid]["family IDs"] |= familyIDs
				seq_table[seqid]["signature_desc"] |= signature_desc

			#this is to assign annotations to GO codes:
			for line in i.split("\n"):
				if '"GO:' in line:
					#this is not stored anywhere?
					GO = re.search(GOpattern, line).group().replace('"GO:', '')
					try:
						signature_desc = re.search(ANNOTpattern, line).group().replace('signature_desc=', '').split(";")[0]
					except AttributeError:
						print("No signature:", line)
						signature_desc = re.search(FAMpattern, line).group().replace('signature_desc=', '').split(";")[0]
					if GO not in GOs_annotation:
						GOs_annotation[GO] = {signature_desc}
					else:
						GOs_annotation[GO].add(signature_desc)

with open("foundGOs.txt", "w") as result:
	print("Writing found GO terms to files...")
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

#reduce reported gene number to gene clusters -> see Jankoviny/fish
with open("IPStable.txt", "w") as result, open("GOtable.txt", "w") as result2,\
open("GO-WEGO.txt", "w") as result3, open("goatools-association.tsv", "w") as result4,\
open("population_from_IPS.tsv", "w") as result5:
	print("Writing IPS table to file...")
	result.write("SeqID\tframe\tRANGE\tGOterms\tIPS ID\tFamily ID\tsignature description\n")
	for k in sorted(seq_table.keys()):
		printframe = seq_table[k]["frame"]
		printrange = seq_table[k]["range"]
		printgo = seq_table[k]["GOs"]
		if len(printgo) == 0:
			printgo = ""
		else:
			gos_goa = ["GO:" + x for x in printgo]
			gos_goa.sort()
			gos = [x for x in printgo]
			gos.sort()
			result2.write("{}\t{}\n".format(k, ",".join(gos)))
			result3.write("{}\t{}\n".format(k, "\t".join(gos_goa)))
			result4.write("{}    {}\n".format(k, ",".join(gos_goa)))
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
		result.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(k, printframe, printrange, printgo, printips, printfamily, printsignature))
		result5.write("{}\n".format(k))

print("parsing finished")