with open("odb9v1_OGs.tsv") as blast_data:
	counter = 0
	for line in blast_data:
		counter += 1
		print(counter + line)
