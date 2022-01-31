import argparse
import os
import gzip
from Bio import SeqIO


def full_clade_set(infilelist):
	"""Prepares a subset of raw (number-only) MATOU seqIDs based on a filtering infile.

	Input:		TSV with seqIDs in the first column

	Returns:	Subset of infile contained in a seqID set

	"""
	clade_set = set()
	for infile in infilelist:
		with open(infile, "rt") as f:
			for x in f:
				clade_set.add(x.split("\t")[0])
				
	return clade_set


def make_clade_set(infilelist):
	"""Prepares a subset of raw (number-only) MATOU seqIDs based on a filtering infile.

	Input:		TSV with seqIDs in the first column

	Returns:	Subset of infile contained in a seqID set

	"""
	clade_set = set()
	for infile in infilelist:
		with open(infile, "rt") as f:
			for x in f:
				seqid = x.split("\t")[0]
				if seqid in pfseqids:
					clade_set.add(seqid)
				
	return clade_set
	

def clades_total_if(dataset="metaT"):
	"""Writes a file containing MATOU abundance per taxon.
	Not to be used, bit complicated to add clades.

	Input:		List of SeqIDs in raw format (number only)

	Returns:	Writes a file

	"""
	all_algae = True
	
	pfseqids = full_clade_set(["Eukaryota.tsv"])

	#sort pfseqids based on taxonomy:	
	pfphaeocystales = full_clade_set(["Phaeocystales.tsv"])
	pfisochrysidales = full_clade_set(["Isochrysidales.tsv"])
	pfhaptophytes = full_clade_set(["Haptophyceae.tsv"])
	if all_algae:
		pfchlorophytes = full_clade_set(["Chlorophyta.tsv"])
		pfcryptophytes = full_clade_set(["Cryptophyceae.tsv"])
		pfdinophytes = full_clade_set(["Dinophyceae.tsv"])
		pfochrophytes = full_clade_set(["Ochrophyta.tsv"])
		pfrhodophytes = full_clade_set(["Rhodophyta.tsv"])
		pfalgae = pfhaptophytes | pfchlorophytes | pfcryptophytes | pfdinophytes | pfochrophytes | pfrhodophytes
	else:
		pfalgae = full_clade_set(["Chlorophyta.tsv", "Cryptophyceae.tsv", "Dinophyceae.tsv", "Haptophyceae.tsv", "Ochrophyta.tsv", "Rhodophyta.tsv"])
	pfother = pfseqids - pfalgae
	
	phaeocystales, isochrysidales, otherhapto, chlorophytes, cryptophytes, dinophytes, ochrophytes, rhodophytes, algae, euks = 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
	total = 0.0
	with gzip.open("{}_occurrences.tsv.gz".format(dataset), "rt") as f,\
		open("{}_{}_if.tsv".format("taxa", dataset), "wt") as result:
		#result.write("PFAM\t#unigenes\tPHC\tISC\tother\ttotal\n")
		if all_algae:
			for l in f:
				seqid = l.split("\t")[0]
				if seqid in pfseqids:
					abundance = float(l.split("\t")[2])
					total += abundance
					if seqid in pfphaeocystales:
						phaeocystales += abundance
					elif seqid in pfisochrysidales:
						isochrysidales += abundance
					elif seqid in pfhaptophytes:
						otherhapto += abundance
					elif seqid in pfchlorophytes:
						chlorophytes += abundance
					elif seqid in pfcryptophytes:
						cryptophytes += abundance
					elif seqid in pfdinophytes:
						dinophytes += abundance
					elif seqid in pfochrophytes:
						ochrophytes += abundance
					elif seqid in pfrhodophytes:
						rhodophytes += abundance
					else:
						euks += abundance
			result.write("Phaeocystales\t{}\nIsochrysidales\t{}\nOtherHaptophyceae\t{}\nChlorophytes\t{}\n".format(phaeocystales, isochrysidales, otherhapto, chlorophytes))
			result.write("Cryptophyceae\t{}\nDinophyceae\t{}\nOchrophyta\t{}\nRhodophyta\t{}\n".format(cryptophytes, dinophytes, ochrophytes, rhodophytes))
			result.write("OtherEuks\t{}\n".format(euks))
			result.write("CtrlTotal\t{}\n".format(total))
		else:
			for l in f:
				seqid = l.split("\t")[0]
				if seqid in pfseqids:
					abundance = float(l.split("\t")[2])
					total += abundance
					if seqid in pfphaeocystales:
						phaeocystales += abund
					elif seqid in phisochrysidales:
						isochrysidales += abund
					elif seqid in pfhaptophytes:
						otherhapto += abund
					elif seqid in pfalgae:
						algae += abund
					else:
						euks += abund
			result.write("Phaeocystales\t{}\nIsochrysidales\t{}\nOtherHaptophyceae\t{}\nOtherAlgae\t{}\n".format(phaeocystales, isochrysidales, otherhapto, algae))
			result.write("OtherEuks\t{}\n".format(euks))
			result.write("CtrlTotal\t{}\n".format(total))
			

def clades_total_dictionary(dataset="metaT"):
	"""Writes a file containing MATOU abundance per taxon.

	Input:		List of SeqIDs in raw format (number only)

	Returns:	Writes a file

	"""
	all_algae = True
	print("Total abundance to be extracted per-clade")
	
	print("  Processing sequence lists...")
	pfseqids = full_clade_set(["Eukaryota.tsv"])
	print("  ...Eukaryota done.")
	#sort pfseqids based on taxonomy:
	pfpelagomonas = full_clade_set(["Pelagomonas.tsv"])
	pfphaeocystales = full_clade_set(["Phaeocystales.tsv"])
	pfisochrysidales = full_clade_set(["Isochrysidales.tsv"])
	pfhaptophytes = full_clade_set(["Haptophyceae.tsv"])
	print("  ...Haptophyceae done.")
	if all_algae:
		pfchlorophytes = full_clade_set(["Chlorophyta.tsv"])
		pfcryptophytes = full_clade_set(["Cryptophyceae.tsv"])
		pfdinophytes = full_clade_set(["Dinophyceae.tsv"])
		pfochrophytes = full_clade_set(["Ochrophyta.tsv"])
		pfrhodophytes = full_clade_set(["Rhodophyta.tsv"])
		pfalgae = pfhaptophytes | pfchlorophytes | pfcryptophytes | pfdinophytes | pfochrophytes | pfrhodophytes
	else:
		pfalgae = full_clade_set(["Chlorophyta.tsv", "Cryptophyceae.tsv", "Dinophyceae.tsv", 
								  "Haptophyceae.tsv", "Ochrophyta.tsv", "Rhodophyta.tsv"])
	print("  ...Algae done.")
	pfother = pfseqids - pfalgae
	
	print("  Reassigning sequence IDs...")
	all_clades = {}
	all_clades.update({x: "OtherEuk" for x in pfother})
	if all_algae:
		all_clades.update({x: "Chlorophyta" for x in pfchlorophytes})
		all_clades.update({x: "Cryptophyceae" for x in pfcryptophytes})
		all_clades.update({x: "Dinophyceae" for x in pfdinophytes})
		all_clades.update({x: "Ochrophyta" for x in pfochrophytes})
		all_clades.update({x: "Rhodophyta" for x in pfrhodophytes})
	else:
		all_clades.update({x: "OtherAlgae" for x in pfalgae})
	#higher to lower rank to reassign
	all_clades.update({x: "OtherHapto" for x in pfhaptophytes})
	#all_clades.update({x: "Isochrysidales" for x in pfisochrysidales})
	#all_clades.update({x: "Phaeocystales" for x in pfphaeocystales})
	all_clades.update({x: "Pelagomonas" for x in pfpelagomonas})
	
	print("  Summing occurrences...")
	abundances = {"Phaeocystales": 0., "Isochrysidales": 0., "OtherHapto": 0., "OtherAlgae": 0., 
	"Chlorophyta": 0., "Cryptophyceae": 0., "Dinophyceae": 0., "Ochrophyta": 0., "Rhodophyta": 0.,  
	"Pelagomonas": 0., 
	"OtherEuk": 0.}
	total = 0.0
	with gzip.open("{}_occurrences.tsv.gz".format(dataset), "rt") as f,\
		open("{}_{}_dict.tsv".format("taxa", dataset), "wt") as result:
		#result.write("PFAM\t#unigenes\tPHC\tISC\tother\ttotal\n")
		for l in f:
			seqid = l.split("\t")[0]
			if seqid in pfseqids:
				abundance = float(l.split("\t")[2])
				total += abundance
				abundances[all_clades[seqid]] += abundance
		if all_algae:
			result.write("Phaeocystales\t{}\nIsochrysidales\t{}\nOtherHaptophyceae\t{}\nChlorophytes\t{}\n".format(
						abundances["Phaeocystales"], abundances["Isochrysidales"], abundances["OtherHapto"], abundances["Chlorophyta"]))
			result.write("Cryptophyceae\t{}\nDinophyceae\t{}\nOchrophyta\t{}\nRhodophyta\t{}\n".format(
						abundances["Cryptophyceae"], abundances["Dinophyceae"], abundances["Ochrophyta"], abundances["Rhodophyta"]))
			result.write("OtherEuks\t{}\n".format(abundances["OtherEuk"]))
			result.write("CtrlTotal\t{}\n".format(total))
		else:				
			result.write("Phaeocystales\t{}\nIsochrysidales\t{}\nOtherHaptophyceae\t{}\nOtherAlgae\t{}\n".format(
						abundances["Phaeocystales"], abundances["Isochrysidales"], abundances["OtherHapto"], abundances["OtherAlgae"]))
			result.write("OtherEuks\t{}\n".format(abundances["OtherEuk"]))
			result.write("CtrlTotal\t{}\n".format(total))


def write_unigenes(infile):
	"""Writes a (nt) subset of the MATOU database that contains unigenes with the given Pfam,
	listed in pfam.tsv.gz. If you intend to translate them for any reason, matching the correct 
	protein translation might be problematic. This is somewhat unnecessary if the database is 
	translated by e.g. Prodigal and searched by HMMER, but comes with other problems (such 
	as setting an appropriate threshold, etc.)

	Input:		List of SeqIDs in raw format (number only)

	Returns:	Writes a file

	"""
	print("Loading list of seqIDs (MATOU format)...")
	with open(infile, "rt") as f:
		fastas = {"MATOU-v1_{}".format(x.split("\t")[0]) for x in f.readlines()}
	print("SeqIDs to keep:", len(fastas))

	print("Searching MATOU fasta db...")
	print("Warning, script will run silently for a long time!")

	suffix = infile.split(".")[-1]
	with open(infile.replace(suffix, "fasta"), "wt") as result,\
		gzip.open("MATOU-v1.fna.gz", "rt") as infasta:
		for seq in SeqIO.parse(infasta, "fasta"):
			if seq.name in fastas:
				result.write(">{}\n{}\n".format(seq.name, seq.seq))


def pfam_occurrences(infile, pfam, dataset="metaT"):
	"""Prepares report files for occurrence_plot.py.

	Input:		Clade-specific TSVs with raw seqIDs in the first column. These were created
				by subsetting the Tara MATOU's taxonomy.tsv.gz using zgrep; modify below to
				work with your target clade.
				Occurrence file provided by Tara MATOU
				A file to extract all seqIDs of a pfam family (FAA-fasta, list of raw seqIDs)
				

	Returns:	File with all occurrences of each seqID across Tara stations (zeroes omitted)
				File with occurrence summaries for each seqID

	"""
	all_algae = True
	
	global pfseqids
	print("Loading list of seqIDs...")
	if infile.split(".")[-1] in ("fasta", "faa", "fas"):
		# Script can process fasta as infile; note that FASTA have modified seqnames! 
		# This is fixed here:
		pfseqids = {x.name.split("_")[1] for x in SeqIO.parse(infile, "fasta")}
	else:
		with open(infile, "rt") as f:
			pfseqids = {x.split("\t")[0] for x in f}

	pfbacteria = make_clade_set(["Bacteria.tsv"])
	pfseqids = pfseqids - pfbacteria

	#sort pfseqids based on taxonomy:	
	pfphaeocystales = make_clade_set(["Phaeocystales.tsv"])
	pfisochrysidales = make_clade_set(["Isochrysidales.tsv"])
	pfhaptophytes = make_clade_set(["Haptophyceae.tsv"])
	if all_algae:
		pfchlorophytes = make_clade_set(["Chlorophyta.tsv"])
		pfcryptophytes = make_clade_set(["Cryptophyceae.tsv"])
		pfdinophytes = make_clade_set(["Dinophyceae.tsv"])
		pfochrophytes = make_clade_set(["Ochrophyta.tsv"])
		pfrhodophytes = make_clade_set(["Rhodophyta.tsv"])
	pfalgae = make_clade_set(["Chlorophyta.tsv", "Cryptophyceae.tsv", "Dinophyceae.tsv", "Haptophyceae.tsv", "Ochrophyta.tsv", "Rhodophyta.tsv"])
	pfother = pfseqids - pfbacteria - pfalgae
	
	all_clades = {}
	all_clades.update({x: "OtherEuk" for x in pfother})
	if all_algae:
		all_clades.update({x: "Chlorophyta" for x in pfchlorophytes})
		all_clades.update({x: "Cryptophyceae" for x in pfcryptophytes})
		all_clades.update({x: "Dinophyceae" for x in pfdinophytes})
		all_clades.update({x: "Ochrophyta" for x in pfochrophytes})
		all_clades.update({x: "Rhodophyta" for x in pfrhodophytes})
	else:
		all_clades.update({x: "OtherAlgae" for x in pfalgae})
	all_clades.update({x: "OtherHapto" for x in pfhaptophytes})
	all_clades.update({x: "Isochrysidales" for x in pfisochrysidales})
	all_clades.update({x: "Phaeocystales" for x in pfphaeocystales})
	
	print("SeqIDS:\n  Phaeocystales: {}\n  Isochrysidales: {}\n  Algae: {}\n  Other: {}".format(len(pfphaeocystales), len(pfisochrysidales), len(pfalgae), len(pfother)))
	print("Warning, script will run silently for a long time!")
	total = 0.0
	phaeocystales_occurrence = 0.0
	isochrysidales_occurrence = 0.0
	other_occurrence = 0.0
	all_occurrences = {}
	with gzip.open("{}_occurrences.tsv.gz".format(dataset), "rt") as f,\
		open("{}_{}-filters.tsv".format(pfam, dataset), "wt") as filters:
		#open("{}_{}.tsv".format(pfam, dataset), "wt") as result:
		#result.write("PFAM\t#unigenes\tPHC\tISC\tother\ttotal\n")
		for l in f:
			try:
				seqid = l.split("\t")[0]
				if seqid in pfseqids:
					total += float(l.split("\t")[2])
					#record the total occurrence of seqid - multiple values are present for more stations
					if seqid not in all_occurrences:
						all_occurrences[seqid] = float(l.split("\t")[2])
					else:
						all_occurrences[seqid] += float(l.split("\t")[2])
					
					filters.write("{}\t{}".format(all_clades.get(seqid, ""), l))
					#lineage-specific occurrences; this is unnecessary as outlier removal occurs later
					#if seqid in pfphaeocystales:
					#	phaeocystales_occurrence += float(l.split("\t")[2])
					#elif seqid in pfisochrysidales:
					#	isochrysidales_occurrence += float(l.split("\t")[2])
					#elif seqid in pfother:
					#	other_occurrence += float(l.split("\t")[2])
			except:
				if l.startswith("unigeneID"):
					pass
				else:
					print("Not a number:", l.strip())
		#result.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pfam, len(pfseqids), phaeocystales_occurrence, isochrysidales_occurrence, other_occurrence, total))

	with open("{}_report.tsv".format(pfam), "wt") as report:
		#report.write("Phaeocystales:\n")
		for seq in pfphaeocystales:
			report.write("Phaeocystales\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
		#report.write("Isochrysidales:\n")
		for seq in pfisochrysidales:
			report.write("Isochrysidales\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
		otherhapto = pfhaptophytes - pfphaeocystales - pfisochrysidales
		for seq in otherhapto:
			report.write("OtherHapto\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
		if all_algae:
			otheralgae = pfalgae - pfhaptophytes - pfchlorophytes - pfcryptophytes - pfdinophytes - pfochrophytes - pfrhodophytes
			for seq in pfchlorophytes:
				report.write("Chlorophyta\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
			for seq in pfcryptophytes:
				report.write("Cryptophyceae\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
			for seq in pfdinophytes:
				report.write("Dinophyceae\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
			for seq in pfochrophytes:
				report.write("Ochrophyta\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
			for seq in pfrhodophytes:
				report.write("Rhodophyta\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))			
		else:	
			otheralgae = pfalgae - pfhaptophytes
		for seq in otheralgae:
			report.write("OtherAlgae\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))
		#report.write("Other:\n")		
		for seq in pfother:
			report.write("OtherEuk\t{}\t{}\n".format(seq, all_occurrences.get(seq, 0)))

	print("Total occurrence: {}".format(total))


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	#use a Pfam ID or a list/fasta of MATOU sequences
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-i', '--infile', help='Pfam file to be used', default='')
	group.add_argument('-p', '--pfam', help='Pfam to analyze', default='')
	parser.add_argument('-c', '--cladestotal', help='Extract clade abundances from MATOU-v1.fna.gz', action='store_true')
	parser.add_argument('-f', '--writefasta', help='Extract listed sequences from MATOU-v1.fna.gz', action='store_true')
	parser.add_argument('-d', '--dataset', help='metaT/metaG', default='metaT')

	args = parser.parse_args()
	if args.cladestotal:
		#this is a Pfam-independent function - collect total occurrence of major taxa.
		clades_total_dictionary(dataset=args.dataset)
	if args.infile != "":
		infile = args.infile
		pfam = infile.split(".")[0]
	elif args.pfam != "":
		infile = args.pfam + ".list"
		pfam = args.pfam
	else:
		#neither pfam nor pfam list was given
		quit("Pfam was not set, finishing!")

	if os.path.isfile(infile) == False:
		os.system('zgrep "{0}" pfam.tsv.gz > {0}.list'.format(pfam))

	if args.writefasta:
		write_unigenes(infile)

	pfam_occurrences(infile, pfam)


if __name__ == '__main__':
	main()

