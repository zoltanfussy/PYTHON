import argparse,os,gzip
from Bio import SeqIO,Entrez
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@natur.cuni.cz'
#update at times:
#ncbi.update_taxonomy_database()


def dictionaries():
	#try using taxids to remove the need to translate the lineage into rank
	global manual_taxids
	manual_taxids = {2461416: ['Bacteria', 'Terrabacteria group', 'Cyanobacteria/Melainabacteria group', 'Cyanobacteria', 'unclassified Cyanobacteria'],
					 3152: ['Eukaryota', 'Viridiplantae', 'Chlorophyta', 'Chloropicophyceae', 'Chloropicales']}

	global eukprot_supergroups
	eukprot_supergroups = {'Alveolata': 'Alveolata', 
						   'Amoebozoa': 'Amoebozoa', 
						   'Ancoracysta': 'basal-euks', 'Ancyromonadida': 'basal-euks', 'Collodictyonidae': 'basal-euks', 'Hemimastigophora': 'basal-euks', 
						   'Malawimonadidae': 'basal-euks', 'Rotosphaerida': 'basal-euks', 
						   'Mantamonas': 'basal-euks', 'Rigifilida': 'basal-euks', 
						   'Telonemia': 'basal-euks',
						   'Apusomonadida': 'basal-opis', 'Apusozoa': 'basal-opis', 'Breviatea': 'basal-opis', 'Filasterea': 'basal-opis', 'Choanoflagellata': 'basal-opis', 
						   'Ichthyosporea': 'basal-opis', 'Pluriformea': 'basal-opis', 'Opisthokonta': 'basal-opis',
						   'Chloroplastida': 'Chloroplastida', 'Viridiplantae': 'Chloroplastida', 
						   'Cryptophyceae': 'Cryptista', 'Kathablepharidacea': 'Cryptista', 'Palpitomonas': 'Cryptista', 
						   'Euglenozoa': 'Discoba', 'Heterolobosea': 'Discoba', 'Jakobida': 'Discoba', 'Tsukubamonadida': 'Discoba',
						   'Fornicata': 'Metamonada', 'Parabasalia': 'Metamonada', 'Preaxostyla': 'Metamonada', 
						   'Fungi': 'Fungi', 
						   'Glaucophyta': 'Glaucophyta', 'Glaucocystophyceae': 'Glaucophyta',
						   'Phaeocystales': 'Haptista-PHC', 'Isochrysidales': 'Haptista-ISO', 'Prymnesiales': 'Haptista-PRY', 'Coccolithales': 'Haptista-COC',
						   'Haptophyta': 'Haptista', 'Centroplasthelida': 'Haptista', 
						   'Metazoa': 'Metazoa', 
						   'Rhodophyta': 'Rhodophyta', 'Picozoa': 'Rhodophyta', 'Rhodelphis': 'Rhodophyta', 
						   'Rhizaria': 'Rhizaria', 
						   'Stramenopiles': 'Stramenopiles',
						   }

	global bacterial_supergroups
	bacterial_supergroups = {'Acidobacteria': 'other-bct', 'Calditrichaeota': 'other-bct', 'Nitrospinae/Tectomicrobia group': 'other-bct', 
							 'Spirochaetes': 'other-bct', 'Synergistetes': 'other-bct', 'Coprothermobacterota': 'other-bct', 
							 'unclassified Bacteria': 'other-bct', 'Bacteria incertae sedis': 'other-bct', 'Nitrospirae': 'other-bct', 
							 'Thermodesulfobacteria': 'other-bct', 'Deferribacteres': 'other-bct', 'Thermotogae': 'other-bct', 'Aquificae': 'other-bct', 
							 'PVC group': 'PVC', 'FCB group': 'FCB-Bacteroidetes', 
							 # subdivisions for Proteobacteria
							 'Alphaproteobacteria': 'Alphaproteobacteria', 'Betaproteobacteria': 'Betaproteobacteria', 'Gammaproteobacteria': 'Gammaproteobacteria',
							 'delta/epsilon subdivisions': 'Deltaproteobacteria',
							 # subdivisions for Terrabacteria
							 'Actinobacteria': 'Actinobacteria', 'Cyanobacteria/Melainabacteria group': 'Cyanobacteria', 'Firmicutes': 'Firmicutes',
							 }

	global viral_supergroups
	viral_supergroups = {'Alphasatellitidae': 'virus', 'Riboviria': 'virus', 'Monodnaviria': 'virus', 'Naldaviricetes': 'virus',
						 #both NCLDV and Polinton viruses belong to Varidnaviria:
						 'Varidnaviria': 'Varidnaviria',
						 'Duplodnaviria': 'Duplodnaviria'
						 }


def read_table(infile):
	outdict = {}
	with open(infile, 'rt') as f:
		for l in f:
			l = l.strip().split("\t")
			if len(l) > 1:
				outdict[l[0]] = int(l[1])
	return outdict


def translate_taxids(taxids_d):
	clade_d = {}
	# 1 = root
	# 131567 = cellular organisms
	# 12908 = unclassified
	# 28384 = other (synthetic, transposons, etc.)
	undefined = {1, 131567, 12908, 28384}
	for key, taxid in taxids_d.items():
		supergroup = ''
		if taxid in undefined:
			continue

		try:
			if taxid in manual_taxids:
				rank = manual_taxids[taxid]
			else:
				lineage = ncbi.get_lineage(taxid)[2:]
				names = ncbi.get_taxid_translator(lineage)
				rank = [names[taxid] for taxid in lineage]
			supergroup = [x for x in rank if x in eukprot_supergroups]
		except ValueError:
			print(f'problem with {key}:{taxid}, please try updating taxa.sqlite')
			with open('errors.log', 'a') as errorfile:
				errorfile.write(f'error retrieving taxid for {key}:{taxid}\n')
			continue

		try:
			# if sequence is bacterial
			if rank[0] == 'Bacteria':
				if len(rank) > 1:
					if rank[1] in ['Terrabacteria group', 'Proteobacteria']:
						if len(rank) > 2:
							supergroup = bacterial_supergroups.get(rank[2], 'other-bct')
						else:
							#print('Insufficient data: {}'.format(rank))
							supergroup = 'other-bct'
					else:
						supergroup = bacterial_supergroups.get(rank[1], 'other-bct')
				else:
					supergroup = 'other-bct'
				#supergroups.append('_'.join(rank[1:2]))

			# if archaeal
			elif rank[0] == 'Archaea':
				supergroup = 'Archaea'

			# if viral
			elif any(x in rank[0] for x in {'viria', 'satellitidae', 'viricetes'}):
				if rank[0] == 'Varidnaviria':
					if rank[2] == 'Nucleocytoviricota':
						supergroup = 'Varidnaviria-NCLDV'
					elif rank[2] == 'Preplasmiviricota':
						supergroup = 'Varidnaviria-Poli'
					else:
						print('Unknown varidnavirus {}'.format(rank))
				else:
					supergroup = viral_supergroups.get(rank[0], 'virus')
				#supergroups.append(rank[0]+'_'+rank[2])
			elif 'viruses' in rank[0]:
				if 'Pleurochrysis sp. endemic' in rank[1]:
					#both Pleurochrysis sp. endemic NCLDV and Polinton virophage
					supergroup = 'Varidnaviria-NCLDV'
				elif 'Pleurochrysis sp. Polinton' in rank[1]:
					supergroup = 'Varidnaviria-Poli'
				else:
					supergroup = 'virus'
				#supergroups.append(rank[0])

			# else eukaryotic
			elif len(supergroup) == 1:
				supergroup = eukprot_supergroups.get(supergroup[0], 'Eukaryota')
			elif len(supergroup) > 1:
				#use the lowest major rank
				supergroup = eukprot_supergroups.get(supergroup[-1], 'Eukaryota')
			#some unwanted situations
			elif rank == ['Eukaryota']:
				supergroup = 'Eukaryota'
			elif supergroup == []:
				print('No usable taxon data: {} => {}'.format(taxid, rank))
				continue
			else:
				print('Other error: {} => {}'.format(taxid, rank))
				continue
			clade_d[key] = supergroup
		except IndexError:
			print("Could not parse taxid", taxid, rank)

	return clade_d


def write_taxids(infile, pfam, taxonomyfile="taxonomy.tsv.gz"):
	"""Compares a subset of the MATOU database that contains unigenes with the given Pfam

	Input:		TSV with seqIDs in the first column

	Returns:	Writes a file

	"""
	print("Loading list of seqIDs...")
	if infile.split(".")[-1] in ("fasta", "faa", "fas"):
		#script can process fasta as infile; note that FASTA have modified seqnames!
		pfseqids = {x.name.split("_")[1] for x in SeqIO.parse(infile, "fasta")}
	else:
		with open(infile, "rt") as f:
			pfseqids = {x.split("\t")[0] for x in f}

	clade_dict = {}
	with gzip.open(taxonomyfile, "rt") as f,\
		open("{}_taxids.tsv".format(pfam), "wt") as taxonfile:
		for l in f:
			seqid = l.split("\t")[0]
			if seqid in pfseqids:
				taxonfile.write("{}\t{}\n".format(seqid, l.split("\t")[1]))


def taxify_file(infile, pfam, translate=True):
	"""Makes dictionaries for translating seqIDs to taxids or clades.

	Input:		TSV with seqIDs in the first column and taxid in the second
				File to modify

	Returns:	Writes a modified file

	"""
	if os.path.isfile("{}_taxids.tsv".format(pfam)) == False:
		print("Prepare *_taxids.tsv first using the -w option!")
	taxid_d = read_table("{}_taxids.tsv".format(pfam))
	print("Imported taxids:", len(taxid_d.keys()))
	suffix = infile.split(".")[-1]
	#Taxids are in raw form - only numbers, therefore in tree and fasta files 
	# we need to extend their names.
	if translate:
		dictionaries()
		outfile = infile.replace(".{}".format(suffix), ".clades.{}".format(suffix))
		clade_d = translate_taxids(taxid_d)
		with open(infile, 'rt') as f:
			file = f.read()
			for seqid in clade_d:
				file = file.replace("MATOU-v1_{}_".format(seqid), "{}.MATOU-v1_{}_".format(clade_d[seqid], seqid))
		with open(outfile, 'wt') as result:
			result.write(file)
		#print("Translated taxids:", clade_d)
	else:
		outfile = infile.replace(".{}".format(suffix), ".taxids.{}".format(suffix))
		with open(infile, 'rt') as f:
			file = f.read()
			for seqid in taxid_d:
				file = file.replace("MATOU-v1_{}_".format(seqid), "{}.MATOU-v1_{}_".format(taxid_d[seqid], seqid))
		with open(outfile, 'wt') as result:
			result.write(file)



def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	#use a Pfam ID or a list/fasta of MATOU sequences
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-i', '--infile', help='Pfam file to be used (fasta or list)', default='')
	group.add_argument('-p', '--pfam', help='Pfam to analyze (file prefix)', default='')
	parser.add_argument('-t', '--taxify', help='Apply taxids in this file (fasta or tree)', default='')
	parser.add_argument('--translate', help='Translate taxids into clades', action='store_true')
	parser.add_argument('-w', '--writetaxids', help='Write a file of taxids', action='store_true')

	args = parser.parse_args()
	if args.infile != "":
		infile = args.infile
		pfam = infile.split(".")[0]
	elif args.pfam != "":
		infile = args.pfam + ".list"
		pfam = args.pfam
	else:
		#neither pfam nor pfam list was given
		quit("Pfam was not set, quitting!")

	#if os.path.isfile(infile) == False:
	#	os.system('zgrep "{0}" pfam.tsv.gz > {0}.list'.format(pfam))

	if args.writetaxids:
		write_taxids(infile, pfam, args.translate)

	taxify_file(args.taxify, pfam, args.translate)

if __name__ == '__main__':
	main()
