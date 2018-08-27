with open("BUSCOlist.txt") as buscolistfile:
	BUSCOlist = set(buscolistfile.read().split("\n"))
outfile = open("busco2og-tsv.txt", "w")
#print(BUSCOlist)

with open("odb9v1_OGs.tsv") as og_data:
	for line in og_data:
		orthogroup = line.split("\t")[0]
		#print(orthogroup)
		if orthogroup in BUSCOlist:
			definition = line.split("\t")[2]
			outfile.write("{}\t{}".format(orthogroup, definition))
