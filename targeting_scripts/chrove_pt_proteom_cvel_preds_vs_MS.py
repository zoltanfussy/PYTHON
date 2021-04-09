from Bio import SeqIO
import re
import os

#THE ONLY VARIABLE IS THE INFILE LIST OF KOs TO PROCESS
home = "/Users/zoliq/ownCloud/"
home = "/Volumes/zoliq data/OwnCloud/"
wd = home + "Terka/Bakalarka/allseq_bestpredictors/sekvence nad-pod threshold/Annotated/"
os.chdir(wd)
with open('pt_proteome_MS.txt') as file:
	MSdata = file.read().split("\n")
	MSdata = set([f.split(".t")[0] for f in MSdata if f != ""])

#fasta = SeqIO.parse("/Users/zoliq/ownCloud/Terka/Bakalarka/allseq_bestpredictors/cvel_FINALSETinclEXT.fasta", "fasta")
fasta = SeqIO.parse(home + "Terka/Bakalarka/allseq_bestpredictors/cvel_FINALSETinclEXT.fasta", "fasta")
fasta_d = {}
seq_d = {}
descpattern = r'gene_product=(.+)\| transcript_product'
for seq in fasta:
	descmatch = re.search(descpattern, seq.description)#.group(1)
	seq_d[seq.name] = seq.seq
	seqname = str(seq.name).split(".t")[0]
	if descmatch:
		desc = descmatch.group(1)
	else:
		desc = "ERROR retrieving annotation"
	if seqname not in fasta_d:
		fasta_d[seqname] = {"desc": desc, "models": [seq.name]}
	elif desc != "ERROR retrieving annotation":
		models = fasta_d[seqname]["models"]
		fasta_d[seqname] = {"desc": desc, "models": models + [seq.name]}

with open("Chromera_signal.txt") as f:
	signal = f.read().split("\n")

with open("cvel_references.txt") as f:
	references = set(f.read().split("\n"))

with open('Chromera_annot_SEKVENCE_ASAFind_nad3.82.txt') as infile:
	cvel_pt_pos_raw = infile.read().split("\n")
	cvel_pt_pos_raw = [f for f in cvel_pt_pos_raw if f != ""]
	cvel_pt_pos = {}
	print("number of predicted positives: {}".format(len(cvel_pt_pos_raw)))
	for i in cvel_pt_pos_raw:
		j = i.split("\t")[0]
		jsplit = j.split(".t")[0]
		if jsplit not in cvel_pt_pos: #zajimaji nas ruzne modely targetovane do stejneho kompartmentu
			cvel_pt_pos[jsplit] = [j]
		else:
			cvel_pt_pos[jsplit].append(j)

with open('Chromera_annot_SEKVENCE_MultiLoc2_nad0.84.txt') as infile:
	cvel_mt_pos_raw = infile.read().split("\n")
	cvel_mt_pos_raw = [f for f in cvel_mt_pos_raw if f != ""]
	cvel_mt_pos = {}
	for i in cvel_mt_pos_raw:
		j = i.split("\t")[0]
		jsplit = j.split(".t")[0]
		if jsplit not in cvel_mt_pos: #zajimaji nas ruzne modely targetovane do stejneho kompartmentu
			cvel_mt_pos[jsplit] = [j]
		else:
			cvel_mt_pos[jsplit].append(j)

with open('Chromera_annot_SEKVENCE_ASAFind_pod3.82.txt') as infile:
	cvel_pt_neg_raw = infile.read().split("\n")
	cvel_pt_neg_raw = [f for f in cvel_pt_neg_raw if f != ""]
	cvel_pt_neg = {}
	for i in cvel_pt_neg_raw:
		j = i.split("\t")[0]
		jsplit = j.split(".t")[0]
		if jsplit not in cvel_pt_neg: #zajimaji nas ruzne modely targetovane do stejneho kompartmentu
			cvel_pt_neg[jsplit] = [j]
		else:
			cvel_pt_neg[jsplit].append(j)

with open('plastNN_cvel.csv') as infile:
	cvel_pt_pos_plastNN_raw = infile.read().split("\n")
	cvel_pt_pos_plastNN_raw = [f for f in cvel_mt_pos_raw if f != ""]
	cvel_pt_pos_plastNN = {}
	for j in cvel_pt_pos_plastNN_raw:
		jsplit = j.split(".t")[0]
		if jsplit not in cvel_pt_pos_plastNN: #zajimaji nas ruzne modely targetovane do stejneho kompartmentu
			cvel_pt_pos_plastNN[jsplit] = [j]
		else:
			cvel_pt_pos_plastNN[jsplit].append(j)

with open('proteom_vs_ASAFind.txt', 'w') as outfile1, \
open('proteom_vs_plastNN.txt', 'w') as outfile2, \
open('plastNN_vs_ASAFind.txt', 'w') as outfile3, \
open('proteom_overlapstat.txt', 'w') as outstats, \
open('proteom_only.fasta', 'w') as outfasta:
	cvel_pt_set = set(cvel_pt_pos)
	#cvel_pt_set = cvel_pt_set - references
	cvel_plastNN_set = set(cvel_pt_pos_plastNN) | references
	print("number of unique positives: ASAFind {}, PlastNN {}".format(len(cvel_pt_set), len(cvel_plastNN_set)))
	print("number of ASAFind vs plastNN overlap: {}".format(len(cvel_pt_set & cvel_plastNN_set)))
	print("number of unique plastNN genes: {}".format(len(cvel_plastNN_set - cvel_pt_set)))
	print("number of unique ASAFind genes: {}".format(len(cvel_pt_set - cvel_plastNN_set)))
	print("_____________________________")
	#print(cvel_pt_set - cvel_plastNN_set)
	#print(cvel_plastNN_set - cvel_pt_set)
	outstats.write("plastNN MS OVERLAP: {}\n".format(len(cvel_plastNN_set & MSdata)))
	outstats.write("plastNN PRED ONLY: {}\n".format(len(cvel_plastNN_set - MSdata)))
	outstats.write("MS minus plastNN: {}\n".format(len(MSdata - cvel_plastNN_set)))
	outstats.write("ASAFind MS OVERLAP: {}\n".format(len(cvel_pt_set & MSdata)))
	outstats.write("ASAFind PRED ONLY: {}\n".format(len(cvel_pt_set - MSdata)))
	outstats.write("MS minus ASAFind: {}\n".format(len(MSdata - cvel_pt_set)))
	outstats.write("plastNN+ASAFind MS OVERLAP: {}\n".format(len((cvel_plastNN_set | cvel_pt_set) & MSdata)))
	outfile1.write("OVERLAP:\n")
	outfile2.write("OVERLAP:\n")

	for i in cvel_pt_set & MSdata:
		outfile1.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
	outfile1.write("\n\n\n\nPRED ONLY:\n")
	for i in cvel_pt_set - MSdata:
		outfile1.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
	outfile1.write("\n\n\n\nMS ONLY:\n")
	for i in MSdata - cvel_pt_set:
		if i in cvel_pt_neg:
			if i in signal:
				local = "signal"
			elif i in cvel_mt_pos:
				local = "mito"
			else:
				local = "cyto"
		else:
			local = "not found"
		try:
			outfile1.write("{}\t{}\t{}\n".format(i, fasta_d[i]["desc"], local))
			for item in fasta_d[i]["models"]:
				outfasta.write(">{}\n{}\n".format(item, seq_d[item]))
		except KeyError:
			outfile1.write("{}\t{}\n".format(i, "pt genome-encoded"))
	
	for i in cvel_plastNN_set & MSdata:
		outfile2.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
	outfile2.write("\n\n\n\nPRED ONLY:\n")
	for i in cvel_plastNN_set - MSdata:
		outfile2.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
	outfile2.write("\n\n\n\nMS ONLY:\n")
	for i in MSdata - cvel_plastNN_set:
		if i in cvel_pt_neg:
			if i in signal:
				local = "signal"
			elif i in cvel_mt_pos:
				local = "mito"
			else:
				local = "cyto"
		else:
			local = "not found"
		try:
			outfile2.write("{}\t{}\t{}\n".format(i, fasta_d[i]["desc"], local))
			#for item in fasta_d[i]["models"]:
			#	outfasta.write(">{}\n{}\n".format(item, seq_d[item]))
		except KeyError:
			outfile2.write("{}\t{}\n".format(i, "pt genome-encoded"))

	for i in cvel_plastNN_set & cvel_pt_set:
		outfile3.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
	outfile3.write("\n\n\n\nplastNN ONLY:\n")
	for i in cvel_plastNN_set - cvel_pt_set:
		#outfile3.write("{}\t{}\n".format(i, fasta_d[i]["desc"]))
		for item in fasta_d[i]["models"]:
			outfile3.write(">{} {}\n \n".format(item, fasta_d[i]["desc"], seq_d[item]))
	outfile3.write("\n\n\n\nASAFind ONLY:\n")
	for i in cvel_pt_set - cvel_plastNN_set:
		for item in fasta_d[i]["models"]:
			outfile3.write(">{} {}\n \n".format(item, fasta_d[i]["desc"], seq_d[item]))

