import os,re,time,argparse
from Bio import SeqIO,Entrez
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@natur.cuni.cz'
#Entrez.api_key = os.environ["API_KEY"]
#print(os.environ)
#update at times:
#ncbi.update_taxonomy_database()

def force_taxid_prot(accession):
	print("Requesting", accession)
	if "|" in accession:
		#this is a bit of a problem
		if len(accession.split("|")) == 4:
			#this happened with old-style headers, now NCBI does not use GI
			accession = accession.split("|")[3]
		else:
			#there should be just one item in accession.split("|") that is not an empty string
			#also, it is supposedly never the first item in such headers
			accession  = [x for x in accession.split("|")[1:] if x not in [""]][0]
		#print(accession)
		
	#Efetch will still fail if a pir| accession is passed!
	try:
		prot = Entrez.efetch(db='protein', id=accession, rettype='gb', retmode='text')
		prot_record = SeqIO.read(prot, 'genbank')
		orgn = prot_record.annotations['organism']
	except:
		print("Could not process", accession)
		taxid = 1
		keeptmpblastfile = True
		with open("errors.log", "a") as errorfile:
			errorfile.write("error retrieving taxid for {}\n".format(accession))
			
	try:
		name2taxid = ncbi.get_name_translator([orgn])
		taxid = name2taxid[orgn][0]
	except:
		print("problem with {}:{}, please try updating taxa.sqlite".format(accession, orgn))
		taxid = 1
		keeptmpblastfile = True
		with open("errors.log", "a") as errorfile:
			errorfile.write("error retrieving taxid for {}:{}\n".format(accession,orgn))

	print("Organism retrieved:", orgn, taxid)

	return taxid

def force_taxid_nucl(accession):
	#lower data traffic with docsum, sufficient for taxid retrieval
	if "|" in accession:
		accession = accession.split("|")[3]
		#print(accession)
	try:
		nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='docsum', retmode='xml')
		record = Entrez.read(nucl)[0]
		taxid = int(record["TaxId"])
		#name2taxid = ncbi.get_name_translator([orgn])
	except:
		print("Could not process", accession)
		keeptmpblastfile = True
		quit()
	print("Organism retrieved:", taxid)

	return taxid

def force_taxid_gbnuc(accession):
	if "|" in accession:
		accession = accession.split("|")[3]
		#print(accession)
	try:
		nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
		nucl_record = SeqIO.read(nucl, 'genbank')
		orgn = nucl_record.annotations['organism']
		name2taxid = ncbi.get_name_translator([orgn])
		taxid = name2taxid[orgn][0]
	except:
		print("Could not process", accession)
		keeptmpblastfile = True
		quit()
	print("Organism retrieved:", orgn, taxid)

	return taxid

############################
###   PARSE PARAMETERS   ###
############################

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print("starting {}".format(current_time))
#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
#dmnd blast input absolutely requires -outfmt 6 qseqid bitscore sseqid qcovhsp pident qlen length !
parser.add_argument('-i', '--infile', help='Diamond outfile(s) to be taxified', required=True)
parser.add_argument('-d', '--work_dir', help='Working directory with the files', default=".")
#parser.add_argument('-a', '--skip_accession_pairing', help='Do not perform parsing prot.accession2taxid', action='store_true')
#input fasta names can be inferred from blast input
parser.add_argument('-nt', '--fasta_nt', help='Nucleotide fasta for parsing', default="")
parser.add_argument('-aa', '--fasta_aa', help='Amino acid fasta for parsing', default="")
parser.add_argument('-g', '--group', help='Taxonomic group to filter', default="Viruses")
parser.add_argument('-n', '--inverse_lookup', help='Inverse lookup', action="store_true")
parser.add_argument('-p', '--pident_threshold', help='percent identity threshold', default=75)
parser.add_argument('-q', '--qcov_threshold', help='query coverage threshold', default=50)
parser.add_argument('-m', '--max_hits', help='Maximum hits to keep', default=100)
parser.add_argument('-t', '--test_mode', help='Testing mode to allow checking tmp files', action="store_true")

args = parser.parse_args()

filetype = "dmnd.out"
transcriptome = "JGI"
if args.infile == "batch": 
	files = [x for x in os.listdir(".") if x.endswith(filetype)]
else:
	files = args.infile.split(",")
print("to analyze:", ", ".join(files))

if args.work_dir != ".":
	os.chdir(args.work_dir)

qthr = float(args.qcov_threshold)
pthr = float(args.pident_threshold)
max_hits = int(args.max_hits)

goodgroups = set(args.group.split(","))
print("Selected taxgroup(s):", goodgroups)
goodgroupsrep = args.group.split(",")[0]

##############
###  MAIN  ###
##############
#quit("Use only for files for which fastas are present in folder")

taxids = {}
try:
	with open("subset.accession2taxid") as taxidfile:
		print("reading accession2taxid...")
		for l in taxidfile:
			l = l.strip().split("\t")#.decode('utf8')
			taxids[l[0]] = l[1]
except:
	print("accession2taxid file not found, skipping")

for i,filepath in enumerate(files):
	path, file = os.path.split(filepath)
	print("\n\n======\nAnalyzing", file)
	# find nucleotide fasta for processing
	if args.fasta_nt == "":
		ntfasta = file.replace(filetype, "fasta")
		if os.path.exists(ntfasta):
			writent = True
		else:
			print("Nucleotide fasta missing/unspecified, nt fasta output muted")
			writent = False
	else:
		ntfasta = args.fasta_nt
		writent = True
	
	# find amino acid fasta for processing
	if args.fasta_aa == "":
		aafasta = file.replace(filetype, "faa")
		if os.path.exists(aafasta):
			writeaa = True
		else:
			print("Amino acid fasta missing/unspecified, aa fasta output muted")
			writeaa = False
			if writent == False:
				print("None of the input fastas provided, skipping!")
				continue
	else:
		aafasta = args.fasta_aa
		writeaa = True
		
	dataset = file.split(".")[0]
	if os.path.isdir("./tmp".format(dataset)) == False:
		os.system("mkdir tmp".format(dataset))
		print("tmp directory created")
	if os.path.isdir("./contaminants".format(dataset)) == False:
		os.system("mkdir contaminants".format(dataset))
		print("contaminant directory created")
	if os.path.isdir("./{}_NR".format(dataset)) == False:
		os.system("mkdir {}_NR".format(dataset))
		print("target directory created")
	statfile = dataset + "_report.txt"
	filt = "tmp/{}_filt.txt".format(dataset)
	check = "tmp/{}_check.txt".format(dataset)
	#use tmpblast file to retrieve taxids if the script crashes on NCBI requests
	#cut -f 3,6 tmp/<blastresult>.tmp >> subset.accession2taxid
	#warning, this file has an atypical column order!
	tmpblast =  "tmp/{}.tmp".format(file)
	keeptmpblastfile = False
	goodscafs = {}
	goodhits = {}
	outscafs = {}
	outhits = {}
	scaf_d = {}
	ranks = {}
	#replacement = {"319938": "288004", "1317118": "1379903", "427920": "1983720"}
	distribution = {goodgroupsrep: 0}
	c = 0 #we need a process monitor
	with open(filepath) as infile, \
	open(tmpblast, "w") as outfile,\
	open(check, "w") as checkfile,\
	open(filt, "w") as filtfile:
		table = infile.read().split("\n")
		#print("{}".format(len(table)))
		print("To be processed: {}".format(len(table)))
		for line in table:
			c += 1
			if c % 10000 == 0:
				print(c)
			if len(line.split("\t")) != 1:
				line = line.split("\t")
				
				#First create a dictionary item
				if len(line) == 8:
					#this is a diamond blastp output, so additional columns can be used for filtering
					score, qcovs, pident = float(line[1]), float(line[3]), float(line[4])
					hitid, sequence = line[2], line[7]
				else:
					print("Not enough columns")
					continue
				try:
					query = line[0]
					#>>>>> MODIFY? <<<<<<#
					query_scaf = query.split("::")[1] # if predicted by transdecoder
					#hitid = line[3]
				except IndexError:
					#assuming protein ID is the same as the contig's! -> prodigal
					query_scaf = query
				if query_scaf not in goodscafs:
					goodscafs[query_scaf] = 0
					outscafs[query_scaf] = 0

				#scaf_d to store best hits
				if query not in scaf_d.keys():
					scaf_d[query] = {"score": score, "qcovs": qcovs, "pident": pident, "taxid": 0, "rank": [""]}
				elif score > scaf_d[query]["score"]:
					#keep the best scoring hit
					scaf_d[query] = {"score": score, "qcovs": qcovs, "pident": pident, "taxid": 0, "rank": [""]}
				#DON'T skip worse hits, unless max_hits is met
				if goodscafs[query_scaf] + outscafs[query_scaf] > max_hits:
					continue

				#Second get taxid and lineage ranks
				if qcovs < 20: #ignore queries not sufficiently covered
					continue
				accession = line[2]
				if accession in taxids:
					taxid = taxids[accession]
				else:
					taxid = force_taxid_prot(accession)
				try:
					lineage = ncbi.get_lineage(taxid)[1:]
					names = ncbi.get_taxid_translator(lineage)
					rank = [names[taxid] for taxid in lineage]
				except ValueError:
					print("problem with {}:{}, please try updating taxa.sqlite".format(accession, taxid))
					keeptmpblastfile = True
					with open("errors.log", "a") as errorfile:
						errorfile.write("error retrieving taxid for {}:{}\n".format(accession, taxid))
					continue

				outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], line[3], line[4], taxid))
				
				if score == scaf_d[query]["score"]:
					scaf_d[query]["taxid"] = taxid
					scaf_d[query]["rank"] = rank
				
				#if any(x in rank for x in goodgroups):
				#	print("Virus!", query)
				if qcovs > qthr and pident > pthr:
					if any(x in rank for x in goodgroups):
						goodscafs[query_scaf] += 1
						goodhits[hitid] = {"taxid": taxid, "sequence": sequence}
					else:
						#higher than the threshold but wrong taxgroup
						outscafs[query_scaf] += 1
						outhits[hitid] = {"taxid": taxid, "sequence": sequence}
				else:
					#lower than the threshold
					outscafs[query_scaf] += 1
					outhits[hitid] = {"taxid": taxid, "sequence": sequence}
				
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

	#print(goodhits)
	#print(goodscafs)
	#now to extract sequence data as requested
	seqlen = 0
	if writent == True:
		in_nucl = SeqIO.parse(ntfasta, "fasta")
		out_nucl = open("{0}_NR/{0}.{1}.NRfilt.fasta".format(dataset, transcriptome), "w")
		out_other = open("contaminants/{0}_other.fasta".format(dataset), "w") #assuming this organism is not in nr
		for seq in in_nucl:
			score = goodscafs.get(seq.name, 0) - outscafs.get(seq.name, 0)
			if args.inverse_lookup == True:
				score = 0 - score
			if score > 0:
				out_nucl.write(">{}\n{}\n".format(seq.name, seq.seq))
				seqlen += len(seq.seq)			
			else:
				#no nr blast hit:
				out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
				seqlen += len(seq.seq)
		print("{} bases extracted".format(seqlen))
		out_nucl.close()
		out_other.close()

	aaseqs = 0
	if writeaa == True:
		in_prot = SeqIO.parse(aafasta, "fasta")
		out_prot = open("{0}_NR/{0}.{1}.NRfilt.faa".format(dataset, transcriptome), "w")
		out_other = open("contaminants/{0}_other.faa".format(dataset), "w")
		for seq in in_prot:
			try:
				query_scaf = seq.name.split("::")[1]
			except IndexError:
				#print("could not parse protID to contig", seq.name)
				query_scaf = seq.name
			score = goodscafs.get(query_scaf, 0) - outscafs.get(query_scaf, 0)
			if args.inverse_lookup == True:
				score = 0 - score

			if score > 0:
				out_prot.write(">{}\n{}\n".format(seq.name, seq.seq))
				aaseqs += 1
			elif seq.name in outscafs:
				out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
			else:
				#ignore:
				pass
				#out_other.write(">{}\n{}\n".format(seq.name, seq.seq))
		out_prot.close()
		out_other.close()
	
	out_hits = open("{0}_NR/{0}.{1}.NRhits.faa".format(dataset, transcriptome), "w")
	for item in goodhits:
		taxid = goodhits[item]["taxid"]
		sequence = goodhits[item]["sequence"]
		out_hits.write(">{}.{}\n{}\n".format(taxid, item, sequence))
	out_hits.close()

	with open(statfile, "w") as result:
		if seqlen != 0:
			result.write("{} bases extracted\n".format(seqlen))
		IL = "NOT " if args.inverse_lookup == True else ""
		print("{}\tsequences of {}{}".format(aaseqs, IL, args.group))
		result.write("{}\tsequences of {}\n".format(aaseqs, goodgroupsrep))
		print("{}\tgood hits".format(len(goodhits.keys())))
		result.write("{}\tgood hits\n".format(len(goodhits.keys())))


	if not args.test_mode:
		os.system("mv {} {}_NR/".format(file, dataset))
		for tmpfile in [filt, check]:
			os.system("rm {}".format(tmpfile))
	if keeptmpblastfile == False:
		os.system("rm {}".format(tmpblast))
	print("now, run quast to analyze contaminants:")
	print("quast {0}_NR/*fasta -o ~/quast/{0}_NR --threads 4".format(dataset))

current_time = time.strftime("%H:%M:%S", t)
print("finished sorting {}".format(current_time))

