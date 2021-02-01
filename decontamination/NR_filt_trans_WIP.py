import os
import re
from Bio import SeqIO,Entrez
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@gmail.com'

def force_taxid(accession):
	#print("Invalid taxid:", accession)
	prot = Entrez.efetch(db='protein', id=accession, rettype='gb', retmode='text')
	prot_record = SeqIO.read(prot, 'genbank')
	orgn = prot_record.annotations['organism']
	name2taxid = ncbi.get_name_translator([orgn])
	try:
		taxid = name2taxid[orgn][0]
	except KeyError:
		taxid = 1
	print("Organism retrieved: {} {}".format(orgn, taxid))

	return taxid

def parse_sseqid(sseqid):
	if sseqid.startswith("gi|"):
		accession = sseqid.split("|")[1]
	else:
		accession = sseqid
	return accession

def get_transcriptID(protID):
	trID = protID.split("::")[1]
	return trID

def counter(c, freq):
	c += 1
	if c % freq == 0:
		print(c)
	return c

##############
###  MAIN  ###
##############
#quit("Use only for files for which fastas are present in folder")
#nepoustet
filetype = "_nr.blastp"
transcriptome = "spades"

files = [x for x in os.listdir(".") if x.endswith(filetype)]
#files = ["bella_tr.Trinity.taxified.out"]
print(files)
goodgroups = {"Metamonada"} #at this point skipping "Preaxostyla" and "Parabasalia" as sister groups, "Excavata" as too general?
goodgroupsrep = "Metamonada"
for file in files:
	print("\n\n======\nAnalyzing", file)
	dataset = file.split(".")[0]
	if not os.path.exists("{0}.{1}.NTfilt.fasta".format(dataset, transcriptome)):
		print("Respective fasta missing!")
		continue
	statfile = dataset + "_report.txt"
	filt = dataset + "_filt.txt"
	check = dataset + "_check.txt"
	cont_bact = set()
	cont_fungal = set()
	cont_animal = set()
	cont_plant = set()
	cont_parabasalia = set()
	cont_other = set()
	species = set()
	goodscafs = {}
	excavs = {}
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
				#First get taxid and lineage rank
				accession = parse_sseqid(line[3])
				trID = get_transcriptID(line[0])
				if line[1] == "N/A":
					taxid = force_taxid(accession)
					lineage = ncbi.get_lineage(taxid)[2:]
				else:
					try:
						taxid = line[1]
						if ";" in line[1]:
							print("Multiple taxids:", line[1])
							taxid = line[1].split(";")[0]
							#for x in line[1].split(";"):
							#	print(ncbi.get_lineage(x)[2:])
						lineage = ncbi.get_lineage(taxid)[2:]
					except ValueError:
						print(f"Invalid taxid, force checking: {accession}")
						taxid = force_taxid(accession)
						lineage = ncbi.get_lineage(taxid)[2:]
				names = ncbi.get_taxid_translator(lineage)
				rank = [names[taxid] for taxid in lineage]
				#Second, prepare ranks to print and make good/bad lists
				if taxid not in species:
					species.add(taxid)
					#orgn = ncbi.get_taxid_translator([taxid])[int(taxid)]
					#print("{}\t{}".format(orgn, "_".join(rank)))
				if any(x in rank for x in goodgroups):
					#print(rank)
					orgn = ncbi.get_taxid_translator([taxid])[int(taxid)]
					goodscafs[trID] = 1
					excavs[trID] = orgn
					distribution[goodgroupsrep] += 1
				elif "Metazoa" in rank:
					#Metazoan ranking is very detailed
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
				else:
					#print(rank)
					group = "_".join(rank[1:3])
					#print(group)
					if group == "":
						#no subgroups defined
						try:
							group = rank[0]
						except IndexError:
							print("WARNING! rank improperly parsed", accession, taxid)
							group = "Bacteria"
					if group not in distribution:
						distribution[group] = 1
					else:
						distribution[group] += 1
				
				#ANY FILTER CAN BE APPLIED
				if len(line) == 6:
					#this is a blastn output, so additional columns can be used for filtering
					bitscore, qcovs, pident = float(line[2]), float(line[4]), float(line[5])
					if qcovs < 50 and pident < 90: #Here i think %id should be higher than nt
						#Sebastian's thresholds:
						#Usually 80% identity and at least 50% coverage of the transcript by the hits. 
						#LGT's don't have such a high identity towards any bacteria usually. 
						goodscafs[trID] = 0
					elif any(x in rank for x in goodgroups):
						goodscafs[trID] = 1
					else:
						#these are high-similarity hits, so sort the sequences as contaminants to subsets:
						if "Bacteria" in rank:
							cont_bact.add(trID)
						elif "Fungi" in rank:
							cont_fungal.add(trID)
						elif "Metazoa" in rank:
							cont_animal.add(trID)
						elif "Streptophyta" in rank:
							cont_plant.add(trID)
						elif "Parabasalia" in rank:
							cont_parabasalia.add(trID)
						else:
							#modify this if organism of interest is in NT
							#if seqID has been added to goodscafs, this will be ignored
							cont_other.add(trID)
							#print(rank)
				else:
					print("Not enough columns")
					continue
				
				#Write ranks to a report
				try:
					checkfile.write(f"{trID}\t{rank[1]}\n")
				except IndexError:
					#this is most likely an unranked bacterium
					#print("\tweird rank", rank)
					checkfile.write(f"{trID}\t{rank}\n")

				#Write ranks of the filtered files
				if trID in goodscafs:
					filtfile.write(f"{trID}\t{rank}\n")


	seqlen = 0
	
	infasta = SeqIO.parse("{0}.{1}.NTfilt.fasta".format(dataset, transcriptome), "fasta")
	outfile = open("{0}_NR/{0}.{1}.NRfilt.fasta".format(dataset, transcriptome), "w")
	#out_bact = open("{0}_NR/{0}_bact.fasta".format(dataset), "w")
	#out_fungal = open("{0}_NR/{0}_fungal.fasta".format(dataset), "w")
	#out_animal = open("{0}_NR/{0}_animal.fasta".format(dataset), "w")
	#out_plant = open("{0}_NR/{0}_plant.fasta".format(dataset), "w")
	#out_para = open("{0}_NR/{0}_para.fasta".format(dataset), "w")
	out_other = open("{0}_NR/{0}_other.fasta".format(dataset), "w") #assuming this organism is not in nr
	for seq in infasta:
		if seq.name in excavs.keys():
			outfile.write(">{} {}\n{}\n".format(seq.name, excavs[seq.name], seq.seq))
			seqlen += len(seq.seq)
		elif seq.name in goodscafs:
			outfile.write(">{}\n{}\n".format(seq.name, seq.seq))
			seqlen += len(seq.seq)
		elif seq.name in cont_bact:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		elif seq.name in cont_fungal:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		elif seq.name in cont_animal:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		elif seq.name in cont_plant:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		elif seq.name in cont_parabasalia:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		elif seq.name in cont_other:
			out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		else:
			#no nr blast hit:
			outfile.write(">{}\n{}\n".format(seq.name, seq.seq))
			seqlen += len(seq.seq)
	print("{} bases extracted".format(seqlen))
	outfile.close()
	#out_bact.close()
	#out_fungal.close()
	#out_animal.close()
	#out_plant.close()
	#out_para.close()
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
	os.system("mv {} {}".format(file, dataset))
	print("now, run quast to analyze contaminants:")
	print("quast {0}_NR/*fasta -o ~/quast/{0}_NR --threads 4".format(dataset))

print("finished sorting")
