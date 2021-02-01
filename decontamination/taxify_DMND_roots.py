import argparse, gzip, time
from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()

Entrez.email = 'zoltan.fussy@gmail.com'

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print("starting {}".format(current_time))
#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Diamond outfile to be taxified', required=True)

args = parser.parse_args()

infile = args.infile

taxids = {}

with open("roots9.accession2taxid") as taxidfile:
	print("reading accession2taxid...")
	for l in taxidfile:
		l = l.strip().split("\t") #.decode('utf8')
		try:
			taxids[l[0]] = l[1]
		except IndexError:
			print(l, "not enough columns")
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)

print("reading diamond output")
with open(infile) as dmndin, open(infile.replace("roots9.dmnd.out", "roots.blastx"), "w") as result:
	for l in dmndin:
		l = l.strip().split("\t")
		#the input format is 
		#qseqid bitscore sseqid qcovhsp pident qlen length
		#default input is:
		#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
		#the output format is:
		#query taxid  bitscore hitID qcovs pident
		sseqid = l[2].split("@")[0]
		taxcode = sseqid.split("_")[0]
		if taxcode in taxids:
			taxid = taxids[taxcode]
		elif sseqid[:8] in taxids:
			taxcode = sseqid[:8] #probably Trichoplax
			taxid = taxids[taxcode]
		else:
			print("code missing", taxcode)
		qcov = l[3] #or alternatively float(l[-1])/float(l[-2])*100
		result.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(l[0], taxid, l[1], sseqid, qcov, l[4]))
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)
print("writing blastx output finished, success!")
