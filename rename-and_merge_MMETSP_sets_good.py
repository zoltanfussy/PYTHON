import os

#os.system("for i in `ls *.tar.gz`;do tar -zxvf $i; done")
renaming = open("renaming_MMETSP.txt").read()
renaming = renaming.split("\n")

codes_d = {}
specieslist = []
skip = ["uncultured eukaryote"]
for line in renaming:
	code = line.split("\t")[0]
	try:
		taxon = line.split("\t")[1]
		specieslist.append(taxon)
	except:
		pass
	codes_d[code] = taxon
multiples = set([])
for species in specieslist:
	if specieslist.count(species) > 1 and species not in skip:
		multiples.add(species)

"""
print(len(multiples))
print(len(set(specieslist)))
"""
unwanted = "doi:10.5281/zenodo.249982-"

filenames = [s for s in os.listdir('.') if s.endswith('.fasta')]

for file in filenames:
	MMETSPcode = file.split(".")[0]
	data = open(file).read().split("\n")
	try:
		species = codes_d[MMETSPcode]
		if species in multiples:
			print("Warning, {} already has a MMETSP code! Merging species datasets...".format(species))
			newfile = ("{}_{}-modif.fasta".format(species, MMETSPcode))
			with open(newfile, 'a') as result:
				for line in data:
					if line.startswith(">"):
						line = line.replace(unwanted,"")
					result.write(line + "\n")
			clustfile = ("{}-clust.fasta".format(species))
			with open(clustfile, 'a') as result:
				for line in data:
					if line.startswith(">"):
						line = line.replace(unwanted,"")
					result.write(line + "\n")
			print("{} added to {}".format(file,clustfile))
		else:
			newfile = ("{}_{}-modif.fasta".format(species, MMETSPcode))
			#os.rename(file, newfile)
			print("{} renamed to {}".format(file,newfile))
			with open(newfile, 'w') as result:
				for line in data:
					if line.startswith(">"):
						line = line.replace(unwanted,"")
					result.write(line + "\n")

	except KeyError as KE:
		print("Not a MMETSP code:", KE)


"""
#more sophisticated
with os.scandir() as it:
	for entry in it:
		if entry.name.startswith('MMETSP') and entry.is_file():
			MMETSPcode = file.split(".")[0]
				try:
					species = codes_d[MMETSPcode]
				except KeyError as KE:
					print("Not a MMETSP code:", KE)
				print(MMETSPcode, species)

unwanted = "doi:10.5281/zenodo.249982-"

filenames = os.listdir('.')
for file in filenames:
	data = open(file).read().split("\n")
	newfile = file.split(".pep")[0] + "-modif.fasta"
	with open(newfile, 'a') as result:
		for line in data:
			if unwanted in line:
				line = line.replace(unwanted,"")
			result.write(line + "\n")
"""
