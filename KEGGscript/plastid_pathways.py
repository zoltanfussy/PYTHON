from Bio import SeqIO
import os
import urllib

home = "/Users/zoliq/ownCloud/"
#home = "/Volumes/zoliq data/OwnCloud/"
wd = home + "Terka/Bakalarka/allseq_bestpredictors"
os.chdir(wd)

orgn = "Vitrella"
compart = "mito" #mito or plastid

# UNCOMMENT TO DOWNLOAD THE CURRENT LIST OF KEGG ENZYMES FROM THE DATABASE
# MIGHT HELP WHEN GETTING KEY ERRORS OF THE EC FILE
"""
url = 'http://rest.kegg.jp/list/ko'
fileName = 'ENZYMES.txt'
urllib.request.urlretrieve(url, fileName)
print("Current list of enzymes downloaded.")
"""
#####################
####     Fx      ####
#####################

def ko_lookup():
	KO_dic = {}
	with open('ENZYMES.txt') as f:
		for line in f:
			ko = line.split()[0][3:]
			#print(ko)
			definition = line.split("\t")[1]
			if '[EC:' in definition:
				EC = line.split('[EC:')[1]
				EC = EC.replace("]","")
				definition = definition.split('[EC:')[0]
				if ' ' in EC:
					EC = EC.split()
				else:
					EC = [EC]
				KO_dic[ko] = [definition, EC]
			else:
				KO_dic[ko] = [definition]
	return KO_dic

def read_pathways(file):
	dic = {}
	order = []
	for line in file:
		if line.startswith("ko"):
			order.append(line.strip())
		elif line[0] in (".", "x", "p", "v"):
			line = line.strip()
			print(line)
			KO = line.split()[1]
			order.append(KO)
			try:
				geneid = line.split()[2].replace(";", "")
			except IndexError:
				geneid = ko_dic[KO][0].split(";")[0]
			try:
				defline = line.split("; ")[1]
			except IndexError:
				defline = ko_dic[KO][0].split("; ")[1]
			if KO not in dic:
				dic[KO] = {"ID": [geneid], "Definition": [defline]}
			else:
				if geneid not in dic[KO]["ID"]:
					dic[KO]["ID"].append(geneid)
				if defline not in dic[KO]["Definition"]:
					dic[KO]["Definition"].append(defline)
	return dic, order

def localize(query):
	with open("sekvence nad-pod threshold/Annotated/Chromera_signal.txt") as f:
		cvelsignal = set(f.read().split("\n"))

	with open("sekvence nad-pod threshold/Annotated/Vitrella_signal.txt") as f:
		vbrasignal = set(f.read().split("\n"))

	with open('sekvence nad-pod threshold/Annotated/Chromera_annot_SEKVENCE_ASAFind_nad3.82.txt') as infile:
		cvel_pt_pos_raw = infile.read().split("\n")
		cvel_pt_pos = {}
		for i in cvel_pt_pos_raw:
			j = i.split("\t")
			if j[0] not in cvel_pt_pos: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				cvel_pt_pos[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Cvel PT positives: {}".format(len(cvel_pt_pos)))

	with open('sekvence nad-pod threshold/Annotated/Chromera_annot_SEKVENCE_MultiLoc2_nad0.84.txt') as infile:
		cvel_mt_pos_raw = infile.read().split("\n")
		cvel_mt_pos = {}
		for i in cvel_mt_pos_raw:
			j = i.split("\t")
			if j[0] not in cvel_mt_pos: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				cvel_mt_pos[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Cvel MT positives: {}".format(len(cvel_mt_pos)))

	with open('sekvence nad-pod threshold/Annotated/Chromera_annot_SEKVENCE_ASAFind_pod3.82.txt') as infile:
		cvel_pt_neg_raw = infile.read().split("\n")
		cvel_pt_neg = {}
		for i in cvel_pt_neg_raw:
			j = i.split("\t")
			if j[0] not in cvel_pt_neg: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				cvel_pt_neg[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Cvel PT negatives: {}".format(len(cvel_pt_neg)))

	with open('sekvence nad-pod threshold/Annotated/Vitrella_annot_SEKVENCE_PredSL_PrediSi_nad0.8.txt') as infile:
		vbra_pt_pos_raw = infile.read().split("\n")
		vbra_pt_pos = {}
		for i in vbra_pt_pos_raw:
			j = i.split("\t")
			if j[0] not in vbra_pt_pos: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				vbra_pt_pos[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Vbra PT positives: {}".format(len(vbra_pt_pos)))

	with open('sekvence nad-pod threshold/Annotated/Vitrella_annot_SEKVENCE_MultiLoc_nad0.88.txt') as infile:
		vbra_mt_pos_raw = infile.read().split("\n")
		vbra_mt_pos = {}
		for i in vbra_mt_pos_raw:
			j = i.split("\t")
			if j[0] not in vbra_mt_pos: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				vbra_mt_pos[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Vbra MT positives: {}".format(len(vbra_mt_pos)))

	with open('sekvence nad-pod threshold/Annotated/Vitrella_annot_SEKVENCE_PredSL_PrediSi_pod0.8.txt') as infile:
		vbra_pt_neg_raw = infile.read().split("\n")
		vbra_pt_neg = {}
		for i in vbra_pt_neg_raw:
			j = i.split("\t")
			if j[0] not in vbra_pt_neg: #nezajimaji nas ruzne modely targetovane do stejneho kompartmentu
				vbra_pt_neg[j[0].split("__")[0]] = (i)
			else:
				pass
	print("Vbra PT negatives: {}".format(len(vbra_pt_neg)))

	outdic = {}

	for item in query:
		if item.startswith("Cvel"):
			if item.split(".t")[0] in cvelsignal:
				SP = " (SP)"
			else:
				SP = ""
			if item in cvel_pt_pos:
				if item in cvel_mt_pos:
					outdic[item] = "pt, also mt" + SP
				elif item in cvel_pt_neg:
					outdic[item] = "pt, also else" + SP
				else:
					outdic[item] = "pt only" + SP
			else:
				if item in cvel_mt_pos:
					outdic[item] = "mt" + SP
				elif item in 	cvel_pt_neg:
					outdic[item] = "else" + SP
				#else:
				#	outdic[item] = "missing",SP
		elif item.startswith("Vbra"):
			if item.split(".t")[0] in cvelsignal:
				SP = " (SP)"
			else:
				SP = ""
			if item in vbra_pt_pos:
				if item in vbra_mt_pos:
					outdic[item] = "pt, also mt" + SP
				elif item in vbra_pt_neg:
					outdic[item] = "pt, also else" + SP
				else:
					outdic[item] = "pt only" + SP
			else:
				if item in vbra_mt_pos:
					outdic[item] = "mt" + SP
				elif item in vbra_pt_neg:
					outdic[item] = "else" + SP
				#else:
				#	outdic[item] = "missing",SP

	return outdic

#####################
####    MAIN     ####
#####################

#prepare KEGG dictionary data
ko_dic = ko_lookup()

#read assumed organellar enzyme list
with open("sekvence nad-pod threshold/Annotated/{}_functions.txt".format(compart)) as f: 
	pathways, order = read_pathways(f)
pwykeys = set(pathways.keys())
print("pathways read")

#read annotations by KAAS
pull_set = set()
#full_set = set() #this is probably not needed
pull_dic = {}
order_data = {}
with open("sekvence nad-pod threshold/Annotated/{}_SEKVENCE2KO.txt".format(orgn)) as f:
	for line in f:
		data = line.strip().split("\t")
		if len(data) == 2:
			#full_set.add(data[0])
			if data[1] in pwykeys:
				pull_set.add(data[0])
				pull_dic[data[0]] = (pathways[data[1]]["ID"], pathways[data[1]]["Definition"], data[1])
				if data[1] not in order_data:
					order_data[data[1]] = [data[0]]
				else:
					order_data[data[1]].append(data[0])
print("list of sequences to extract prepared")

preds = localize(pull_set)

#filter fastas
badchars = ("*?UOJBZX")
if orgn == "Chromera":
	infasta = SeqIO.parse("cvel_FINALSETinclEXT.fasta", "fasta")
elif orgn == "Vitrella":
	infasta = SeqIO.parse("vbra_FINALSETinclEXT.fasta", "fasta")
with open("sekvence nad-pod threshold/Annotated/{}_{}set.fasta".format(orgn, compart), "w") as out:
	for seq in infasta:
		if seq.name in pull_set:
			sequence = "".join(a for a in str(seq.seq) if a not in badchars)
			if "X" in sequence:
				sequence.replace("X", "")
				print("error", seq.name)
			out.write(">{} {}; {}\n{}\n".format(seq.name, pull_dic[seq.name][0][0], " ".join(pull_dic[seq.name][1]), sequence))

print("sequences parsed")

#now write results into a table
with open("sekvence nad-pod threshold/Annotated/{}_{}set.tsv".format(orgn, compart), "w") as out2:
	for i in order:
		if i.startswith("ko"):
			out2.write("\n{}\n===\n".format(i))
		else:
			for s in order_data.get(i, ["NONE"]):
				pred = preds.get(s.split("__")[0], "NONE")
				out2.write("{}\t{}; {}\t{}\t{}\n".format(i, pathways[i]["ID"][0], " ".join(pathways[i]["Definition"]), s, pred))

print("outtable written, quitting")
