import gzip

clades = set()
with gzip.open("taxonomy.tsv.gz", "rt") as handle:
	for l in handle.read():
		l = l.split("\t")
		if len(l) == 4:
			clades.add(l[2])

with open("taxonomies.txt", "wt") as result:
	result.write("\n".join(list(clades)))
