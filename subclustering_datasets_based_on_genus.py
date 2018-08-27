import os

filenames = [f for f in os.listdir('.') if os.path.isfile(f)]
mergedfile = '.merge.fas'

genera = {}
skip = ["unwanted_species"]

for file in filenames:
	if file.endswith("fasta"):
		genus = file.split("_")[0]
		if genus in genera:
			genera[genus].append(file)
		else:
			genera[genus] = [file]
counter = 0
for genus in genera:
	if len(genera[genus]) > 1 and genus not in skip:
		newfile = genus + mergedfile
		for file in genera[genus]:
			towrite = open(newfile, 'a')
			towrite.write(open(file).read())
		print(newfile + " being clustered")
		counter+=1
		os.system("usearch -cluster_fast {} -id 0.95 -centroids {}.clus.fasta -notrunclabels".format(newfile, genus))
	else:
		newfile = genera[genus][0]
		print(newfile + " being clustered")
		counter+=1
		os.system("usearch -cluster_fast {} -id 0.95 -centroids {}.clus.fasta -notrunclabels".format(newfile, genus))

os.system("rm *.merge.fas")
print(counter, "files written")
os.system("for i in `ls *clus.fasta`; do mv $i clustered; done")