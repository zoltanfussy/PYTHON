import argparse, gzip, time
from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()

Entrez.email = 'zoltan.fussy@gmail.com'

def force_taxid(accession):
	print("WARNING: missing taxid in input file, requesting from NCBI server:", accession)
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

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print("starting {}".format(current_time))
#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Diamond outfile to be taxified', required=True)

args = parser.parse_args()

infile = args.infile

accessionlist = set()
with open(infile) as dmndin:
	c = 0
	for l in dmndin:
		l = l.strip().split("\t")
		#the input format is 
		#qseqid bitscore sseqid qcovhsp pident qlen length
		accessionlist.add(l[2])
		c += 1
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)
print("parsing blastp output finished")
print("{} unique records found".format(c))

taxids = {}
with gzip.open("prot.accession2taxid.gz", mode='rb') as taxidfile, open("subset.accession2taxid", "w") as subset:
#with open("subset.accession2taxid") as taxidfile, open("subset", "w") as subset:
	print("reading accession2taxid...")
	c = 0
	for l in taxidfile:
		#c += 1
		#if c % 10000000 == 0:
		#	print("{}M".format(c/1000000))
		l = l.strip().split("\t")#.decode('utf8')
		try:
			if l[1] in accessionlist:
				subset.write("{}\n".format("\t".join(l)))
				taxids[l[1]] = l[2]
		except IndexError:
			print(l, "not enough columns")
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)
print("reading accession2taxid finished")

with open(infile) as dmndin, open(infile.replace("out", "blastp"), "w") as result:
	for l in dmndin:
		l = l.strip().split("\t")
		#the input format is 
		#qseqid bitscore sseqid qcovhsp pident qlen length
		#default input is:
		#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
		#the output format is:
		#query taxid  bitscore hitID qcovs pident
		sseqid = l[2]
		if sseqid in taxids:
			taxid = taxids[sseqid]
		elif sseqid == "*":
			taxid = 1
		else:
			taxid = force_taxid(sseqid)
		qcov = l[3] #or alternatively float(l[-1])/float(l[-2])*100
		result.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(l[0], taxid, l[1], sseqid, qcov, l[4]))
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)
print("writing blastp output finished, success!")
