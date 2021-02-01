from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()

"""with open("fetch_entrez_missing.txt") as f:
	missing = f.read().split("\n")
#print(missing)
name2taxid = ncbi.get_name_translator(missing)
print(name2taxid)
#TESTING PURPOSES ONLY:
#missing = ['Homo', 'Aspergillus', 'Haloquadratum']
#taxid2name = ncbi.get_taxid_translator([9606, 9443])
#rankofnode = ncbi.get_rank([9606, 9443])

for genus in name2taxid:
	#retriev at least one species:
	descendants = ncbi.get_descendant_taxa(genus)
	lineage = ncbi.get_lineage(descendants[0])[2:7]
	names = ncbi.get_taxid_translator(lineage)
	rank = [names[taxid] for taxid in lineage]
	if "Eukaryota" in rank:
		rank.remove("Eukaryota")
	print("{}\t{}\t{}".format(genus, ";".join(rank)))"""
data_d = {}
with open("/Users/morpholino/OwnCloud/genomes/EukProt_included_data_sets.v02.2020_06_30.txt") as f:
#with open("/Users/morpholino/OwnCloud/AndyLab/MMETSP_mcl/species_list.txt") as f:
	missing = f.read().split("\n")
	if missing[0].startswith("EukProt_ID"):
		missing = missing[1:]
	for x in missing:
		data_d[x.split("\t")[1].replace("_", " ")] = {"EukProt_ID": x.split("\t")[0], "Data": "\t".join(x.split("\t")[2:])}
	missing = [x.split("\t")[1].replace("_", " ") for x in missing] #this contains species taxids now
print(missing)
name2taxid = ncbi.get_name_translator(missing)
print(name2taxid)

for species in missing: #name2taxid:
	if species in name2taxid:
		continue
	#print(species)
	#retriev at least one species:
	try:
		genus = species.split()[0]
		descendants = ncbi.get_descendant_taxa(genus)
		lineage = ncbi.get_lineage(descendants[0])[2:] #[2:7]
		names = ncbi.get_taxid_translator(lineage)
		for taxid in lineage:
			if names[taxid] == genus:
				#print(species, taxid)
				name2taxid[species] = [taxid]
		rank = [names[taxid] for taxid in lineage]
	except:
		print(species, "not found")
		name2taxid[species] = ["N/A"]
	#if "Eukaryota" in rank:
	#	rank.remove("Eukaryota")
	#print("{}\t{}".format(genus, ";".join(rank))) #name2taxid[genus][0], 

with open("/Users/morpholino/OwnCloud/genomes/EukProt_included_data_sets.v02.taxids.txt", "w") as result:
	for item in name2taxid:
		data = data_d[item]
		result.write("{}\t{}\t{}\t{}\n".format(data["EukProt_ID"], name2taxid[item][0], item.replace(" ", "_"), data["Data"]))


