import os

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "progs/SCRIPTS/Kinannote_1.0/results/"
os.chdir(wd)

filelist = [f for f in os.listdir(".") if f.endswith("txt")]
print("Files to be analyzed: ", filelist, "\n\n")

#first, create a list of all kinases in the results
nonredundant = set()

for item in filelist:
	with open(item) as file:
		kinases = file.read().split('\n')
	if kinases[-1] == '':
		kinases = kinases[:-1]
	for kinase in kinases:
		#print(kinase)
		name = kinase.split('\t')[0]
		nonredundant.add(name)

nonredundant = list(nonredundant)
nonredundant.sort()
print("List of found kinases: ", nonredundant, "\n\n\n==============================")

#here analyze data and populate table
kinase_d = {}
n = 0

for item in filelist:
	if n == 0:
		kinase_d["taxon"] = [item.split('.txt')[0]]
	else:	
		kinase_d["taxon"].append(item.split('.txt')[0])
	with open(item) as file:
		kinases = file.read().split('\n')
		if kinases[-1] == '':
			kinases = kinases[:-1]
		#print(kinases)
		item_kinases = {}
		for kinase in kinases:
			#print(kinase)
			name = kinase.split('\t')[0]
			count = kinase.split('\t')[1]
			item_kinases[name] = count
		print(item, item_kinases.keys())
		for kinase in nonredundant:
			#print(kinase)
			if kinase in item_kinases.keys():
				try:
					kinase_d[kinase].append(item_kinases[kinase])
				except KeyError:
					kinase_d[kinase] = [item_kinases[kinase]]
			else:
				try:
					kinase_d[kinase].append('')
				except KeyError:
					kinase_d[kinase] = ['']
	n += 1

#print(kinase_d)
print("\n\n====================\nTable of kinases prepared, now writing...")

with open("results_table.tsv", "w") as tsv_file:
	for key in kinase_d.keys():
		counts = '\t'.join(kinase_d[key])
		tsv_file.write("{}\t{}\n".format(key, counts))

print("Writing file finished.")
