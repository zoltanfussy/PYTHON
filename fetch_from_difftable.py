import os
from Bio import SeqIO,Entrez
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
Entrez.email = 'zoltan.fussy@gmail.com'

#ncbi.update_taxonomy_database()

def read_blast_results(filetype):
	seq_d = {}
	fastas = [x for x in os.listdir(".") if x.endswith("fasta")]
	for file in fastas:
		for seq in SeqIO.parse(file, "fasta"):
			seqname = parse_sseqid(seq.name)
			seq_d[seqname] = seq.seq

	blast_d = {}
	files = [x for x in os.listdir(".") if x.endswith(filetype)]
	for file in files:
		with open(file) as f:
			for l in f:
				l = l.strip().split("\t")
				query = l[0]
				target = parse_sseqid(l[1])
				try:
					blast_d[query] = {}
					blast_d[query]["target"] = target
					blast_d[query]["sequence"] = seq_d[target]
				except KeyError:
					#fasta sequence missing
					pass

	return blast_d

def get_seq_chromosome(accession):
	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
	nucl_record = SeqIO.read(nucl, 'fasta')
	acc = nucl_record.name
	sequence = nucl_record.seq
	print("SEQ retrieved:", acc)

	return sequence

def get_aa(accession):
	handle = Entrez.efetch(db='protein', id=accession, rettype='fasta', retmode='text')
	record = SeqIO.read(handle, 'fasta') #only one record is expected
	#nucl_sequence = record['GBSeq_sequence']
	acc = record.name
	sequence = record.seq
	print("SEQ retrieved:", acc)

	return sequence

def get_aa_from_nucl(accession):
	#explanation of db structure:
	"""nucl = Entrez.einfo(db='nucleotide')
	record = Entrez.read(nucl)
	record.keys()
	record["DbInfo"]"""

	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='xml')
	record = Entrez.read(nucl)[0] #only one record is expected
	#nucl_sequence = record['GBSeq_sequence']
	features = record['GBSeq_feature-table']
	features = [x for x in features if x['GBFeature_key'] == 'CDS'][0] #only one CDS will be extracted
	#orgn = nucl_record.annotations['organism']
	for item in features['GBFeature_quals']:
		if item['GBQualifier_name'] == 'protein_id':
			prot_id = item['GBQualifier_value']
		if item['GBQualifier_name'] == 'translation':
			prot_sequence = item['GBQualifier_value']
			break
	print("SEQ retrieved:", accession)
	sequence = prot_id, prot_sequence #if CDS should be translated

	return sequence

def get_aa_from_gene(accession):
	#to retrieve the protein id, the easiest is to parse the gene table
	handle = Entrez.efetch(db='gene', id=accession, rettype='gene_table', retmode='text')
	gene_table = handle.read()
	gene_table = gene_table.split("\n")
	prot_id = ""
	for line in gene_table:
		if line.startswith("protein"):
			#maybe use re?
			prot_id = [x for x in line.split() if x.startswith("XP")][0].replace(",", "")
			#prot_id = line.split()[1].replace(",", "")
			print("ACC retrieved:", prot_id)
			if not prot_id.startswith("XP"):
				print(line)
			break

	if prot_id != "":
		handle = Entrez.efetch(db='protein', id=prot_id, rettype='fasta', retmode='text')
		record = SeqIO.read(handle, 'fasta')
		sequence = record.seq
		print("SEQ retrieved:", prot_id)
	else:
		sequence = ""

	#for xml - more data available:
	#explanation of db structure:
	"""nucl = Entrez.einfo(db='gene')
	record = Entrez.read(nucl)
	record.keys()
	record["DbInfo"]"""
	#handle = Entrez.efetch(db='gene', id=accession, rettype='gb', retmode='xml')
	#record = Entrez.read(handle)[0]
	#print(record.keys())
	#sequence accessions are in ['Entrezgene_comments']


	return prot_id, sequence

def parse_sseqid(sseqid):
	if sseqid.startswith("gi|") or sseqid.startswith("gb|"):
		accession = sseqid.split("|")[1]
	else:
		accession = sseqid
	return accession

##############
###  MAIN  ###
##############
#quit("Use only for files for which fastas are present in folder")
#nepoustet
filetype = "met-gtf.tsv"
files = [x for x in os.listdir(".") if x.endswith(filetype)]
#files = ["bella_tr.Trinity.taxified.out"]
print(files)

blast_d = read_blast_results("blastn")
#print(blast_d)

#try accessions:
#XM_019096574.1	> refseq nucl column 15
#LOC109045000	> gene locus column 12
#NC_031698.1	12720705	12720796	>chromosome coordinates columns 0,1,2

for file in files:
	print("\n\n======\nAnalyzing", file)
	accessionset = set()
	protset = set()
	processedfile = file.replace(filetype, "ok")
	if os.path.exists(processedfile):
		print("record file exists")
		with open(processedfile) as f:
			for l in f:
				l = l.strip().split()
				accessionset.add(l[3])
				if l[3].startswith("XP"):
					protset.add("{}_{}-{}".format(l[0], l[1], l[2]))
			print(len(accessionset), "items finished")
	dataset = file.split(".")[0]
	processed = open(processedfile, "a")
	aa_result = open(file.replace(filetype, "faa"), "a")
	nt_result = open(file.replace(filetype, "fna"), "a")
	with open(file) as f:
		for l in f:
			line = l.strip().split("\t")
			if line[-1] == "protein_id":
				continue
			query = "{}_{}-{}".format(line[0], line[1], line[2]) #chromosome
			if query in protset: #we already have a protein for this locus, skip
				continue
			if line[17].startswith("XP"):
				accession = line[17]
				if accession in accessionset:
					continue
				sequence = get_aa(accession)
				accessionset.add(accession) #only after the fetch, accession is added
				protset.add(accession)
				print("{}\t{}\t{}\t{}\trecorded".format(line[0], line[1], line[2], accession))
				processed.write("{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], accession))
				aa_result.write(">{} query={}:{}-{}\n{}\n".format(accession, line[0], line[1], line[2], sequence))
			elif line[12].startswith("LOC"):
				accession = line[12] #will be 12
				if accession in accessionset:
					continue
				aa_accession, sequence = get_aa_from_gene(accession)
				accessionset.add(accession)
				protset.add(accession)
				print("{}\t{}\t{}\t{}\trecorded".format(line[0], line[1], line[2], accession))
				processed.write("{}\t{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], accession, aa_accession))
				aa_result.write(">{} query={}:{}-{}\n{}\n".format(aa_accession, line[0], line[1], line[2], sequence))
			elif line[15].startswith("XM"):
				accession = line[15] #will be 15
				if accession in accessionset:
					continue
				aa_accession, sequence = get_aa_from_nucl(accession) # protein AC will not be retrieved - fix issue
				accessionset.add(accession)
				protset.add(accession)
				print("{}\t{}\t{}\t{}\trecorded".format(line[0], line[1], line[2], accession))
				processed.write("{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], accession, aa_accession))
				aa_result.write(">{} query={}:{}-{}\n{}\n".format(aa_accession, line[0], line[1], line[2], sequence))
			else:
				start, end = int(line[1]), int(line[2])
				if query in accessionset:
					continue
				if query in blast_d:
					sequence = blast_d[query]["sequence"]
					accession = blast_d[query]["target"]
					print("SEQ stored:", accession)
				else:
					sequence = get_seq_chromosome(line[0])
					sequence = sequence[start+1:end]
					accession = query
				accessionset.add(query)
				print("{}\t{}\t{}\t{}\trecorded as NT".format(line[0], line[1], line[2], accession))
				processed.write("{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], query))
				nt_result.write(">{} query={}\n{}\n".format(accession, query, sequence))

			if len(sequence) > 10000:
				print("long transcript!", accession)
	processed.close()
	aa_result.close()
	nt_result.close()
			
print("Finished")

