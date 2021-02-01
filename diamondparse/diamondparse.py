import os, re, argparse

###################
# parse arguments #
###################
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--pool', help='Pool results?', default=True)
parser.add_argument('-q', '--query', help='Blast query file', default='diamondquery.fasta')
parser.add_argument('-r', '--results', help='Comma separated blast result file(s)', default='diamond.out')

args = parser.parse_args()
outpool = args.query.split(".")[0]
results = args.results.split(",") #results is always a list!
print(results)

def seqparse(fastafile):
	seqdict = {}
	with open(fastafile) as f:
		data = f.read()
		data = data.split(">")[1:]
		for seq in data:
			seqname = seq.split("\n")[0]
			seqname = seqname.split(" ")[0]
			seqseq = "".join(seq.split("\n")[1:])
			seqdict[seqname] = {seqname: seqseq}
	return seqdict

###################
#  go to workdir  #
###################

if not os.path.isdir("diamond_out"):
	os.system("mkdir diamond_out")
	print("DIAMONDPARSE: Creating target directory")

print("DIAMONDPARSE: Starting...")
#FILENAMES:
SORTDIR = 'diamond_out/'
QUERYFILE = args.query
LOG = 'Diamond-mod-log.txt' #NOTE this is for a different than default-setting diamond outfiles
#evaluethresh = 0.0001 #e-value not in output

###################
#      MAIN       #
###################

#create seq dictionary based on query file
sequence_dict = seqparse(QUERYFILE)
seq_count = len(list(sequence_dict.keys()))
print("DIAMONDPARSE: Sequences extracted from query file: {}".format(seq_count))

#goodhits = set()
#badqueries = set()

#search for taxon in brackets
taxonpattern = r'\[(.+)\]'
for BLASTOUT in results:
	with open(BLASTOUT) as infile:
		for line in infile:
			line = line.strip().split("\t")
			query = line[0]
			hitname = line[1]
			hitdesc = line[-2]
			if hitname[:6] in ("contig", "isotig"):
				#this happens with some of the MMETSP sequences
				if line[-2] == "Phaeocystis antarctica, Strain CCMP1374":
					hitname = "Phaant-CCMP1374_" + hitname
				elif hitname.endswith("ChaeIllumina"):
					hitname = "Chaetoceros-UNC1202_" + hitname.replace("-ChaeIllumina", "")
				elif hitname.endswith("NitzIllumina"):
					hitname = "Nitzschia-sp.-RCC80_" + hitname.replace("-NitzIllumina", "")
			if BLASTOUT == "nr.out":
				#modify fasta header if description available
				#NOTE: for REFSEQ_BLASTOUT there are no descriptions
				if hitname != hitdesc:
					taxon = re.search(taxonpattern, hitdesc).group(1)
					if taxon:
						hitdesc = hitdesc.replace("PREDICTED: ", "")
						taxon = "_".join(taxon.split()[:2])
						hitdesc = hitdesc.split(" [")[0]
						hitname = "{}_{}".format(taxon, hitdesc)
			if query not in sequence_dict:
				print("ERROR! unspecified query:", query)
			else:
				sequence_dict[query].update({hitname: line[-1]})

if args.pool in {True, "True", "true"}:
	print("DIAMONDPARSE: pooling hit sequences to {}{}...".format(SORTDIR, outpool))
	uniqseqs = set()
	with open("{}{}_out.fasta".format(SORTDIR, outpool), "w") as result: 
		for query in sequence_dict:
			for item in sequence_dict[query]:
				#print(item)
				hitseq = sequence_dict[query][item]
				#print(hitseq)
				#print(item, hitseq)
				if hitseq not in uniqseqs:
					uniqseqs.add(hitseq)
					result.write(">{}\n{}\n".format(item, hitseq))

else:
	print("DIAMONDPARSE: pooling disabled, hit sequences in {}".format(SORTDIR))
	for query in sequence_dict:
		uniqseqs = set()
		with open("{}{}_out.fasta".format(SORTDIR, query.replace("/", "_")), "w") as result: 
			for item in sequence_dict[query]:
				#print(item)
				hitseq = sequence_dict[query][item]
				#print(hitseq)
				#print(item, hitseq)
				if hitseq not in uniqseqs:
					uniqseqs.add(hitseq)
					result.write(">{}\n{}\n".format(item, hitseq))

print("DIAMONDPARSE: Extracted proteins written to files. Finished.")
