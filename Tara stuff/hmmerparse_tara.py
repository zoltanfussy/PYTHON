import os,argparse
from Bio import SeqIO
import gzip
import itertools

#os.chdir("/Users/zoliq/ownCloud/progs/PYTHON/mydata")
"""
prefixes = [x for x in os.listdir("./alignments") if x.endswith(".ali")]
#prefixes = ["PF00134_rp75", "PF02984_rp75"]
for prefix in prefixes:
	#os.system("hmmbuild {0}.hmm {0}.ali".format(prefix))
	os.system("hmmsearch -o {0}.result.txt --notextw --tblout {0}.table.txt --cpu 2 {0}.hmm hmmdb.fa".format(prefix, prefix))
"""

def make_clade_set(infilelist):
	"""Compares a subset of the MATOU database that contains unigenes with the given Pfam

	Input:		TSV with seqIDs (number-only) in the first column

	Returns:	Subset of infile contained in a seqID set

	"""
	if isinstance(infilelist, str):
		infilelist = infilelist.split(",")
	clade_set = set()
	for infile in infilelist:
		with open(infile, "rt") as f:
			for x in f:
				seqid = x.split("\t")[0]
				if seqid in hmmids_nt:
					clade_set.add(seqid)
				
	return clade_set
	

# def hmm_hits(genefile):
# 	#provide a set of all hmmer hits
# 	if genefile == "batch":
# 		files = [x for x in os.listdir(".") if x.endswith("table.txt")]
# 	else:
# 		files = genefile.split(",")
# 	files.sort()
# 	hmmids_nt = set()
# 	for file in files:
# 		with open(file) as f:
# 			for line in f:
# 				if "inclusion threshold" in line or "Domain annotation" in line:
# 					print("HMMER inclusion threshold reached, wrapping up...")
# 					break
# 				line = line.split()
# 				if len(line) < 9:
# 					continue
# 				if  line[0] != "#":
# 					item = line[8]
# 					if not item.startswith("MATOU"):
# 						print("Problem with line", " ".join(line))
# 					#modify item so that prodigal appendix is removed ==> gene ID
# 					#item = item[:item.rfind("_")]
# 					hmmids_nt.add(item)
# 
# 	return hmmids_nt


def hmm_hits_thresholded(genefile):
	"""Provides a set of thresholded hmmer hits; 
	Expects one file from function "occurrence" but can work with multiple, 
	e.g. to process entire pathways.

	Input		List of files as a ","-separated string

	Returns 	Set of seqIDs from HMMER hmmids_nt_filtered.

	"""
	if genefile == "batch":
		files = [x for x in os.listdir(".") if x.endswith("table.txt")]
	else:
		files = genefile.split(",")
	files.sort()
	multidomain = set()
	written = {} #move multidomain, written into the for loop to have domain-specific files
	for file in files:
		with open(file) as infile:
			for line in infile:
				if "inclusion threshold" in line or "Domain annotation" in line:
					print("HMMER inclusion threshold reached, wrapping up...")
					break
				line = line.split()
				if len(line) < 9:
					continue
				if  line[0] != "#":
					try:
						evalue = float(line[0])
						#seqname = line[8][:line[8].rfind("_")] #trim the prodigal appendix
						seqname = line[8].split("_")[1] #extract raw seqID
						#print(seqname)
						if seqname in written:				
							pass
						elif evalue < 1: #which is the usual case for hits above inclusion threshold
							written[seqname] = line[8]
							multidomain.add(seqname)
						elif seqname in multidomain: #another weak domain found in the same seq
							written[seqname] = line[8]
						else: #when there is a weak domain in the target
							multidomain.add(seqname)
					except ValueError:
						pass
						#this is not a score table line
	return written


def writeseqs(fastadb, genefile, lineagefile):
	"""Writes a subset of MATOU database that were found as HMM hits

	Input:		List of SeqIDs

	Returns:	Writes a file

	"""
	#prepare a file list
	if genefile == "batch":
		files = [x for x in os.listdir(".") if x.endswith("table.txt")]
	else:
		files = genefile.split(",")
	files.sort()
	print("Files to be processed", ", ".join(files))

	global hmmids_nt
	hmmids_nt = hmm_hits_thresholded(genefile)
	if lineagefile != "":
		lineage_subset = make_clade_set(lineagefile)
		if len(lineage_subset) > 0:
			hmmids_nt_filtered = { x: hmmids_nt[x] for x in hmmids_nt.keys() if x in lineage_subset }
	else:
		hmmids_nt_filtered = hmmids_nt
	hmmids_aa = set(hmmids_nt_filtered.values())
	print(type(hmmids_aa))
		
	print("Total of {} hits found".format(len(hmmids_nt)))
	if lineagefile != "":
		print("Filtered hits: {}".format(len(hmmids_nt_filtered)))
	print("Examples:", [x for x in itertools.islice(hmmids_nt_filtered, 10)])
	print("Example seqids:", [x for x in itertools.islice(hmmids_aa, 10)])

	#filter seqs identified as hits
	print("Now to extract sequences from db file", fastadb)
	if len(files) == 1:
		#single file mode
		file = files[0]
		print("Processing for FASTA extraction:", file)
		outname = file.replace(".table.txt", "-dbhits.fasta")
		with open(file) as infile, open(outname, "w") as result:
			if fastadb.split(".")[-1] in ("fasta", "fas", "faa"):
				fastafile = SeqIO.parse(fastadb, "fasta")
				for seq in fastafile:
					if seq.name in hmmids_aa:
						result.write(">{}\n{}\n".format(seq.name, str(seq.seq).replace("*", "")))
			elif fastadb.endswith("gz"):
				with gzip.open(fastadb, "rt") as unzipped:
					for seq in SeqIO.parse(unzipped, "fasta"):
						if seq.name in hmmids_aa:
							result.write(">{}\n{}\n".format(seq.name, str(seq.seq).replace("*", "")))
			else:
				quit("fasta database invalid format")

	else:
		seq_d = {}
		print("Warning, script will run silently for a long time!")
		if fastadb.split(".")[-1] in ("fasta", "fas", "faa"):
			fastafile = SeqIO.parse(fastadb, "fasta")
			for seq in fastafile:
				if seq.name in hmmids_aa:
					seq_d[seq.name] = str(seq.seq).replace("*", "")
		elif fastadb.endswith("gz"):
			with gzip.open(fastadb, "rt") as fastafile:
				for seq in SeqIO.parse(fastafile, "fasta"):
					if seq.name in hmmids_aa:
						seq_d[seq.name] = str(seq.seq).replace("*", "")
		else:
			quit("fasta database invalid format")
		print("Total of {} sequences matched hits".format(len(seq_d.keys())))	

		#filter results and export fasta
		for file in files:
			print("Processing for FASTA extraction:", file)
			multidomain = set()
			written = set() #move multidomain, written into the for loop to have domain-specific files
			name = file.replace(".table.txt", "-dbhits.fasta")
			with open(file) as infile, open(name, "w") as result:
				for line in infile:
					if "inclusion threshold" in line or "Domain annotation" in line:
						print("HMMER inclusion threshold reached, wrapping up...")
						break
					line = line.split()
					if len(line) < 9:
						continue
					if  line[0] != "#":
						#seqnames match between prodigal faa and hmmer hit hmmids_aa
						seqname = line[8]
						if seqname in hmmids_aa:				
							result.write(">{}\n{}\n".format(seqname, seq_d[seqname]))
	print("Sequence extraction finished!")


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-f', '--infile', help='Fasta database file to process', default='/mnt/mokosz/home/zoli/proj/MATOU/MATOU-v1.prodigal.faa.gz')
	parser.add_argument('-l', '--lineage', help='Lineage filter', default='')
	parser.add_argument('-d', '--directory', help='Workdirectory', default='.')
	parser.add_argument('-g', '--genefile', help='Gene file', default='batch')
	parser.add_argument('-e', '--export_hits', help='Write HMMER hits as fasta', action='store_true')

	args = parser.parse_args()
	if args.directory != ".":
		os.chdir(args.directory)
	if args.export_hits:
		writeseqs(args.infile, args.genefile, args.lineage)


if __name__ == '__main__':
	main()

