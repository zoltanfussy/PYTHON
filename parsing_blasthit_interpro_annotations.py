import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_protein
import argparse
import re

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

def predictseq(positive, ipshitrange, ipsannotation):
	nonsensemutation = False
	# checking if query is more than 10 aa shorter than best hit, flagging accordingly
	j = ipshitrange[0]
	k = ipshitrange[1]
	fromto = [j,k]
	orfN = "full"

	# finding encoded protein, first upstream
	stops = ["*"]
	aalist = []
	downstream = allseq_d[positive]
	aacount = 0
	#when in upstream region, sequences are reset when stop codon is met
	while aacount < j:
		aa = downstream[:1]
		aacount += 1
		downstream = downstream[1:]
		aalist.append(aa)
		if aa in stops:
			aa = ''
			aalist = []
	
	#M is assumed as start codon in upstream region
	if "M" in aalist:
		start = aalist.index("M")
		protein_seq = sum(aalist[start:], Seq("", generic_protein))
		newfrom = fromto[0] - len(protein_seq) + 1
	else:
		protein_seq = sum(aalist, Seq("", generic_protein))
		newfrom = fromto[0] - len(protein_seq) + 1

	# finding encoded protein, now entering conserved region with the IPS hit
	minimumquerylen = len(protein_seq) + k
	while len(protein_seq) < minimumquerylen and len(downstream) > 0:
		aa = downstream[:1]
		if aa in stops:
			protein_seq = protein_seq + "X"
			nonsensemutation = True
		else:
			protein_seq = protein_seq + aa
		downstream = downstream[1:]

	# finding encoded protein, downstream of the homology region
	while aa not in stops and len(downstream) > 0:
		aa = downstream[:1]
		protein_seq = protein_seq + aa
		downstream = downstream[1:]

	#print(protein_seq[-3:])
	if aa in stops:
		orfC = "full"
	else:
		orfC = "partial"
	fromto = (newfrom, newfrom + len(protein_seq))
	#print("{} processed".format(sequence))

	if nonsensemutation:
		flag = "in-frame non-sense codon"
	elif orfC == "full" and orfN == "full":
		flag = "complete"
	elif orfC == "full" and orfN == "partial":
		flag = "N-partial"
	elif orfC == "partial" and orfN == "full":
		flag = "C-partial"
	else:
		flag = "fragment"

	fullseq = allseq_d[positive]
	ipsid = ipsannotation.split(":")[0]
	newitem = {"Bestorg": "IPShit", "Besthit": ipsid, "Sequence": protein_seq, 
			"Flag": flag, "Range": fromto, "IPSannot": ipsannotation}
	return newitem

##################################
#### Create working directory ####
##################################

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
if True:
	print("changing to default directory")
	defdir = "Jankoviny/Tick_transcriptome/"
	wd = home + defdir
	os.chdir(wd)


#FILENAMES:
QUERYFILE = SeqIO.parse('Trinity_all_trimmed_stages-AA.fasta', "fasta")
FILTHITS = 'NR_filtered_hits.txt'
FILTANNOTS = 'NR_filtered_hits_annotated.txt'
FILTFASTA = SeqIO.parse('NR_goodproteins.fasta', "fasta")
BLASTHEADERS = "blastheaders.txt"
IPSTABLE = "IPStable.txt"

#parse BLAST headers first
print("Parsing nr descriptions into annotations")
#accessionpattern = r'[_A-Z0-9]{6,12})\.1'
#look up pattern with:
#ACCID = re.search(accessionpattern, description).group(1)

superfluous = {"AltName:", "Flags:", "TPA:", " ", "", "X2", "X3", "X4"}
blastannot_d = {}
c = 0
with open(BLASTHEADERS) as f:
	for l in f:
		c += 1
		descriptions = l.strip().replace(">", "").split("] ")
		taxa = set()
		keywords = []
		key = descriptions[0].split(" ")[0]
		for item in descriptions:
			try:
				orgn = item.replace("]", "").split("[")[1]
				description = item.split("[")[0]
				for word in description.split(" ")[1:]: #first "word" is the accession number
					if word in superfluous:
						break
					word = word.replace("Full=", "").replace("Short=", "").replace(";", "")
					if word not in keywords:
						keywords.append(word)
				taxa.add(orgn)
			except IndexError:
				for word in item.split(" ")[1:]: #first "word" is the accession number
					if word in superfluous:
						break
					word = word.replace("Full=", "").replace("Short=", "").replace(";", "")
					if word not in keywords:
						keywords.append(word)
				#print(item.split()[0], "item has no orgn but was added to annotations")

		blastannot_d[key] = {"organisms": ",".join(list(taxa)), "descriptions": " ".join(keywords)}
		#FOR TESTING ONLY
		#if key == "5ELS_A":
		#	print("this should be stored")
		#	print(blastannot_d[key])
		#if c == 1000:
		#	quit("demo finished")


#link FILTHITS with their annotations
print("Linking filtered hits with their annotations")
with open(FILTHITS) as f:
	data = f.read()
entries = data.split("==================\n")
entries.pop() #the last entry is empty
blasthits = set()
fasta_db = {}
#this can be done only once:
"""
print("Writing updated table of hits with annotations")
with open(FILTANNOTS, "w") as result, open("no_annot_fetch.txt", "w") as error:
	for entry in entries:
		query = entry.split("\n")[0]
		if query != "":
			query = query.split()[0]
			#print(query)
			blasthits.add(query)
		hits = entry.split("\n")[2:]

		for item in hits:
			try:
				hitacc = item.split()[0]
				#FOR TESTING ONLY:
				#if hitacc == "5ELS_A":
				#	print("this should be loaded")
				#	print(blastannot_d[hitacc])
				if hitacc[-2] == ".":
					hitacc = hitacc.split(".")[0] + ".1"

				description = blastannot_d.get(hitacc, "no description!")
				#print(description)
				if description == "no description!":
					print(hitacc, "no description ")
					print(entry)
					error.write("{}\n".format(hitacc))
					result.write("{}\tDESC={} ORGN={}\n".format(item, "no description retrieved", "NONE"))
				else:
					result.write("{}\tDESC={} ORGN={}\n".format(item, description["descriptions"], description["organisms"]))
			except IndexError:
				pass
quit()
"""
print("Looking up and merging annotations from all hits per query")
c = 0
total = len(entries)

print(total, "queries to process")
for entry in entries:
	c += 1
	if c % 5000 == 0:
		print("{:.2f}%".format(c/total*100))
	query = entry.split("\n")[0]
	if query != "":
		query = query.split()[0]
		#print(query)
		blasthits.add(query)
	hits = entry.split("\n")[2:]
	tophitid = hits[0].split()[0]
	if tophitid[-2] == ".":
		tophitid = tophitid.split(".")[0] + ".1"
	try:
		tophitorg = blastannot_d[tophitid]["organisms"]
	except KeyError:
		tophitorg = "no description"
	#print(tophitid, tophitorg)
	keywordspool = []
	for item in [hits[0]]: #change to "hits" to proceed through all hits
		try:
			hitacc = item.split()[0]
			#FOR TESTING ONLY:
			#if hitacc == "5ELS_A":
			#	print("this should be loaded")
			#	print(blastannot_d[hitacc])
			if hitacc[-2] == ".":
				hitacc = hitacc.split(".")[0] + ".1"

			description = blastannot_d.get(hitacc, "no description!")
			#print(hitacc, description)
			if description == "no description!":
				if keywordspool == []: 
					keywordspool = ["no description retrieved"]
				#print(hitacc, "no description ")
				#result.write("{}\tDESC={} ORGN={}\n".format(item, "no description retrieved", "NONE"))
				#pass
			else:
				d = description["descriptions"].split(" ")
				for word in d:
					if word not in keywordspool:
						keywordspool.append(word)
				#result.write("{}\tDESC={} ORGN={}\n".format(item, description["descriptions"], description["organisms"]))
				#pass
		except IndexError:
			pass
	if "no description retrieved" in keywordspool:
		if len(keywordspool) > 1:
			keywordspool.remove("no description retrieved")
	if "UNKNOWN" in keywordspool:
		if len(keywordspool) > 5:
			keywordspool.remove("UNKNOWN")
	fasta_db[query] = {"Besthit": tophitid, "Bestorg": tophitorg, "Descriptions": " ".join(keywordspool), "Range": ""}
print("100.00%")

print("Extracting seq ranges, flags")
for seq in FILTFASTA:
	rangedata = seq.description.split("RANGE=")[1]
	seqflag = seq.description.split(" RANGE=")[0].split("FLAG=")[1]
	query = seq.name
	seqrange = (int(rangedata.split("-")[0]),int(rangedata.split("-")[1]))
	fasta_db[query]["Range"] = seqrange
	fasta_db[query]["Sequence"] = seq.seq
	fasta_db[query]["Flag"] = seqflag

print("Hit assignment finished")
print("Loading original query sequences for IPS ORF extraction")
allseq_d = {}
for seq in QUERYFILE:
	allseq_d[seq.name] = seq.seq

#find if seqs have overlapping IPS / BLAST annotations:
print("Finding overlaps between IPS and BLAST annotations")
notin = 0
with open(IPSTABLE) as f:
	for l in f:
		if not l.startswith("SeqID"):
			l = l.strip()		
			contig = l.split()[0]
			hr = l.split("\t")[2].split("-")
			ipshitrange = (int(hr[0]),int(hr[1]))
			ipsannotation = "{}:{}".format(l.split("\t")[5], l.split("\t")[6])
			if len(contig.split("_")) == 4:
				positive = contig
			else:
				positive = contig + "_" + l.split()[1]
			#print(positive)
			if positive in blasthits:
				blasthitrange = fasta_db[positive]["Range"]
				#print(positive, ipshitrange, fasta_db[positive])
				if overlaps(ipshitrange, blasthitrange) in ("embed", "overlap"):
					fasta_db[positive].update({"IPSannot": ipsannotation})
				else: 
					print("non-overlapping annotations: {} {} vs {}".format(positive, ipshitrange, blasthitrange))
					notin += 1
					newitem = predictseq(positive, ipshitrange, ipsannotation)
					#print(positive, newitem)
					fasta_db[positive + "_IPS"] = newitem
			else:
				notin += 1
				print("new sequence annotation", positive)
				newitem = predictseq(positive, ipshitrange, ipsannotation)
				#print(positive, newitem)
				fasta_db[positive] = newitem

print("{} IPS-annotated sequences not among nr hits!".format(notin))
print("Writing final annotated fasta file")
writefastas = list(fasta_db.keys())
writefastas.sort()
#print(writefastas)
total = len(writefastas)
c = 0
with open("NR_goodproteins_final.fasta", "w") as writefile:
	for item in writefastas:
		c += 1
		if c % 5000 == 0:
			print("{:.2f}%".format(c/total*100))
		writedata = fasta_db[item]
		#print(item, writedata)
		if writedata.get("Descriptions", "") != "":
			annotations = writedata["Descriptions"]
		else:
			annotations = ""
		if writedata.get("IPSannot", "") != "":
			annotations = "{} {}".format(annotations, writedata["IPSannot"])
		#add flag here, adjust range format
		seqdesc = "DESC={} FLAG={} RANGE={}-{} BBH={} ORG={}".format(annotations, writedata["Flag"], writedata["Range"][0], writedata["Range"][1], writedata["Besthit"], writedata["Bestorg"])
		writefile.write(">{} {}\n{}\n".format(item, seqdesc, writedata["Sequence"]))
print("100.00%")
print("Parsing finished")