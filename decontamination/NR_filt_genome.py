import os,re
from Bio import SeqIO,Entrez
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@natur.cuni.cz'
#Entrez.api_key = os.environ["API_KEY"]
#update at times:
#ncbi.update_taxonomy_database()

def force_taxid_prot(accession):
	print("Requesting", accession)
	if "|" in accession:
		accession = accession.split("|")[3]
		#print(accession)
	prot = Entrez.efetch(db='protein', id=accession, rettype='gb', retmode='text')
	prot_record = SeqIO.read(prot, 'genbank')
	orgn = prot_record.annotations['organism']
	name2taxid = ncbi.get_name_translator([orgn])
	taxid = name2taxid[orgn][0]
	print("Organism retrieved:", orgn, taxid)

	return taxid

def force_taxid_nucl(accession):
	if "|" in accession:
		accession = accession.split("|")[3]
		#print(accession)
	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
	nucl_record = SeqIO.read(nucl, 'genbank')
	orgn = nucl_record.annotations['organism']
	name2taxid = ncbi.get_name_translator([orgn])
	taxid = name2taxid[orgn][0]
	print("Organism retrieved:", orgn, taxid)

	return taxid

def parse_sseqid(sseqid):
	if sseqid.startswith("gi|"):
		accession = sseqid.split("|")[1]
	else:
		accession = sseqid
	return accession

def del_pX(nodename):
	partpattern = r'p\d+'
	try:
		partID = re.search(partpattern, nodename).group()
		nodename = nodename.replace(partID, "")
		#this also works
		#for hit in re.findall(partpattern, nodename):
		#	nodename = nodename.replace(hit, "")
	except:
		pass
		#print(f"No pattern in {nodename}")
	return nodename

def counter(c, freq):
	c += 1
	if c % freq == 0:
		print(c)
	return c

#previously fetched from NCBI
lookups = {1267768: [2, 1224, 28211, 204455, 31989, 1855413, 1267768], 
		1379903: [2, 1224, 28211, 204455, 31989, 93682, 1379903], 
		93064: [2, 1224, 28211, 204457, 41297, 13687, 93064], 
		1514140: [2759, 33634, 2696291, 2836, 589449, 33837, 232508, 210587, 1529381, 1514140], 
		2494563: [28883, 10699, 196894, 2494563], 
		1333996: [2, 1224, 28211, 356, 41294, 1649510, 1333996], 
		301: [2, 1224, 1236, 72274, 135621, 286, 136841, 1232139, 301],
		2605946: [2, 1224, 28211, 204455, 31989, 58840, 2605946]}

##############
###  MAIN  ###
##############
#quit("Use only for files for which fastas are present in folder")
#nepoustet
filetype = "dmnd.blastx"
transcriptome = "spades"

files = [x for x in os.listdir(".") if x.endswith(filetype)]
#files = ["bella_tr.Trinity.taxified.out"]
print(files)
goodgroups = {"Stramenopila", "Stramenopiles"} #at this point skipping "Preaxostyla" and "Parabasalia" as sister groups, "Excavata" as too general?
goodgroupsrep = "Stramenopila"
for file in files:
	print("\n\n======\nAnalyzing", file)
	if not os.path.exists(file.replace(filetype, "fasta")):
		print("Respective fasta missing!")
		continue
	dataset = file.split(".")[0]
	if os.path.isdir("./{}_NR".format(dataset)) == False:
		os.system("mkdir {}_NR".format(dataset))
		print("target directory created")
	statfile = dataset + "_report.txt"
	filt = dataset + "_filt.txt"
	check = dataset + "_check.txt"
	cont_bact = list()
	cont_fungal = list()
	cont_animal = list()
	cont_plant = list()
	cont_parabasalia = list()
	cont_other = list()
	species = set()
	goodscafs = {}
	outscafs = {}
	scaf_d = {}
	ranks = {}
	#replacement = {"319938": "288004", "1317118": "1379903", "427920": "1983720"}
	distribution = {goodgroupsrep: 0}
	c = 0 #we need a process monitor
	with open(file) as infile, \
	open(check, "w") as checkfile,\
	open(filt, "w") as filtfile:
		table = infile.read().split("\n")
		#print("{}".format(len(table)))
		print("To be analyzed: {}".format(len(table)))
		for line in table:
			c += 1
			if c % 10000 == 0:
				print(c)
			if len(line.split("\t")) != 1:
				line = line.split("\t")
				
				#First get taxid and lineage ranks
				accession = parse_sseqid(line[3])
				taxid = line[1]
				if taxid == "N/A":
					taxid = force_taxid_prot(accession)
					lineage = ncbi.get_lineage(taxid)[2:]
				elif taxid in lookups:
					lineage = lookups[taxid]
				else:
					try:
						if ";" in taxid:
							print("Multiple taxids:", taxid)
							taxid = taxid.split(";")[0]
							#for x in taxid.split(";"):
							#	print(ncbi.get_lineage(x)[2:])
						lineage = ncbi.get_lineage(taxid)[2:]
					except ValueError:
						print(f"Invalid taxid, force checking: {line} {accession}")
						taxid = force_taxid_prot(accession)
						lineage = ncbi.get_lineage(taxid)[2:]
						lookups[taxid] = lineage
				names = ncbi.get_taxid_translator(lineage)
				rank = [names[taxid] for taxid in lineage]
				
				#Second, create a dictionary item
				try:
					query = line[0]
					query_scaf = del_pX(query)
					if query_scaf not in goodscafs:
						goodscafs[query_scaf] = 0
						outscafs[query_scaf] = 0
					score = float(line[2])
					#hitid = line[3]
				except IndexError:
					#this should not happen, but testing what other errors might happen
					continue

				if query not in scaf_d.keys():
					scaf_d[query] = {"taxid": taxid, "score": score}
				elif score > scaf_d[query]["score"]:
					#keep the best scoring hit
					scaf_d[query] = {"taxid": taxid, "score": score}
				else:
					continue
				
				#Third, prepare ranks to print and store distribution data
				if taxid not in species:
					species.add(taxid)
					#orgn = ncbi.get_taxid_translator([taxid])[int(taxid)]
					#print("{}\t{}".format(orgn, "_".join(rank)))
				if any(x in rank for x in goodgroups):
					#print(rank)
					orgn = ncbi.get_taxid_translator([taxid])[int(taxid)]
					ranks[query_scaf] = orgn
					group = goodgroupsrep
					distribution[goodgroupsrep] += 1
				elif "Metazoa" in rank:
					#Metazoan ranking is too detailed
					try:
						#print(rank[7])
						group = "Opisthokonta_Metazoa_" + rank[7]
						if group not in distribution:
							distribution[group] = 1
						else:
							distribution[group] += 1
					except IndexError:
						if not "Trichoplax" in rank:
							print("out of range:" + "_".join(rank))
					ranks[query_scaf] = group
				else:
					#print(rank)
					group = "_".join(rank[1:3])
					#print(group)
					if group == "":
						#no subgroups defined
						try:
							group = rank[0]
						except IndexError:
							group = "N/A"
					if group not in distribution:
						distribution[group] = 1
					else:
						distribution[group] += 1
					ranks[query_scaf] = group
				
				#ANY FILTER CAN BE APPLIED
				if len(line) == 6:
					#this is a blastn output, so additional columns can be used for filtering
					bitscore, qcovs, pident = float(line[2]), float(line[4]), float(line[5])
					if qcovs < 50 and pident < 80: #Sebastian's thresholds:
						#Usually 80% identity and at least 50% coverage of the transcript by the hits. 
						#LGT's don't have such a high identity towards any bacteria usually. 
						goodscafs[query_scaf] += 1
						if "Bacteria" in rank:
							cont_bact.append(query_scaf)
							outscafs[query_scaf] += 1
					elif any(x in rank for x in goodgroups):
						goodscafs[query_scaf] += 1
					else:
						#these are high-similarity hits, so sort the sequences as contaminants to subsets:
						if "Bacteria" in rank:
							cont_bact.append(query_scaf)
							outscafs[query_scaf] += 1
						elif "Fungi" in rank:
							cont_fungal.append(query_scaf)
							outscafs[query_scaf] += 1
						elif "Metazoa" in rank:
							cont_animal.append(query_scaf)
							outscafs[query_scaf] += 1
						elif "Streptophyta" in rank:
							cont_plant.append(query_scaf)
							outscafs[query_scaf] += 1
						elif "Parabasalia" in rank:
							cont_parabasalia.append(query_scaf)
							outscafs[query_scaf] += 1
						else:
							#modify this if organism of interest is in NT
							#if seqID has been added to goodscafs, this will be ignored
							cont_other.append(query_scaf)
							outscafs[query_scaf] += 1
							#print(rank)
				else:
					print("Not enough columns")
					continue
				
				#Write ranks to a report
				try:
					checkfile.write(f"{query}\t{rank[1]}\n")
				except IndexError:
					#this is most likely an unranked bacterium
					#print("\tweird rank", rank)
					checkfile.write(f"{query}\t{rank}\n")

				#Write ranks of the filtered files
				if query in goodscafs:
					filtfile.write(f"{query}\t{rank}\n")

	seqlen = 0
	
	infasta = SeqIO.parse(file.replace(filetype, "fasta"), "fasta")
	outfile = open("{0}_NR/{0}.{1}.NRfilt.fasta".format(dataset, transcriptome), "w")
	out_bact = open("{0}_NR/{0}_bact.fasta".format(dataset), "w")
	out_fungal = open("{0}_NR/{0}_fungal.fasta".format(dataset), "w")
	out_animal = open("{0}_NR/{0}_animal.fasta".format(dataset), "w")
	out_plant = open("{0}_NR/{0}_plant.fasta".format(dataset), "w")
	out_para = open("{0}_NR/{0}_para.fasta".format(dataset), "w")
	out_other = open("{0}_NR/{0}_other.fasta".format(dataset), "w") #assuming this organism is not in nr
	for seq in infasta:
		if goodscafs.get(seq.name, 0) > outscafs.get(seq.name, 0):
			outfile.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
			seqlen += len(seq.seq)			
			#out_para.write("{}:\t{}\t{}\n".format(seq.name, goodscafs.get(seq.name, 0), outscafs.get(seq.name, 0)))
		#there might be a lot of low-coverage hits, but most are from bacteria
		elif cont_bact.count(seq.name) > goodscafs.get(seq.name, 0):
			out_bact.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		#note that the following are only high-similarity, high-coverage hits
		elif cont_fungal.count(seq.name) > goodscafs.get(seq.name, 0):
			out_fungal.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		elif cont_animal.count(seq.name) > goodscafs.get(seq.name, 0):
			out_animal.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		elif cont_plant.count(seq.name) > goodscafs.get(seq.name, 0):
			out_plant.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		elif cont_parabasalia.count(seq.name) > goodscafs.get(seq.name, 0):
			out_para.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		elif cont_other.count(seq.name) > goodscafs.get(seq.name, 0):
			out_other.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
		else:
			#no nr blast hit:
			outfile.write(">{} {}\n{}\n".format(seq.name, ranks.get(seq.name, ""), seq.seq))
			seqlen += len(seq.seq)
			#out_para.write("{}:\t{}\t{}\n".format(seq.name, goodscafs.get(seq.name, 0), outscafs.get(seq.name, 0)))
	print("{} bases extracted".format(seqlen))
	outfile.close()
	out_bact.close()
	out_fungal.close()
	out_animal.close()
	out_plant.close()
	out_para.close()
	out_other.close()
	
	print("{} different species as hits".format(len(species)))
	with open(statfile, "w") as result:
		if seqlen != 0:
			result.write("{} bases extracted\n".format(seqlen))
		groups = list(distribution.keys())
		groups.sort()

		g = goodgroupsrep
		print("{}\tsequences of {}".format(distribution[g], g))
		result.write("{}\tsequences of {}\n".format(distribution[g], g))
		groups.remove(goodgroupsrep)

		for g in groups:
			print("{}\tsequences of {}".format(distribution[g], g))
			result.write("{}\tsequences of {}\n".format(distribution[g], g))
	os.system("mv {} {}_NR/".format(file, dataset))
	print("now, run quast to analyze contaminants:")
	print("quast {0}_NR/*fasta -o ~/quast/{0}_NR --threads 4".format(dataset))

print("NCBI taxid lookups", lookups)
print("finished sorting")
