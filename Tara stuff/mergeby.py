def _subset(string, option):
	result = ""
	if string[2].isnumeric():
		options = {"station": slice(0,3),
				   "layer": slice(3,6),
				   "filter": slice(7,11)}
	elif string[1].isnumeric():
		options = {"station": slice(0,2),
				   "layer": slice(2,5),
				   "filter": slice(6,10)}
		#add zero if two-digit station number
		if "station" in option:
			result = "0"		
	else:
		options = {"station": slice(0,1),
				   "layer": slice(1,4),
				   "filter": slice(5,9)}
		#add two zeroes if one-digit station number
		if "station" in option:
			result = "00"

	for o in option:
		result += string[options[o]]

	return result


def mergeby(file, outfile, by="all"):
	if by == "all":
		options = ["station", "layer", "filter"]
	else:
		by = by.split(",")
		#station always first
		options = [x for x in ["station", "layer", "filter"] if x in by]
		if any(x not in ("station", "filter", "layer") for x in by):
			unknown = [x for x in by if x not in ("all", "station", "filter", "layer")]
			print("! unrecognized merging option omitted: {}".format(unknown))


	codes = {'CCKK': '0.22-3µm', 
			 'CCII': '0.2-1.6µm', 
			 'GGMM': '0.8-5µm', 
			 'KKQQ': '3-20µm', 
			 'MMQQ': '5-20µm', 
			 'GGZZ': '0.8-2000µm', 
			 'QQSS': '20-180µm', 
			 'QQRR': '20-200µm',
			 'SSUU': '180-2000µm', 
			 }

	merged = {}

	if file.endswith("gz"):
		import gzip
		with gzip.open(file, "rt") as f:
			for l in f:
				l = l.split()
				if len(l) == 4:
					subset = _subset(l[2], options)
					occur = float(l[3])
					if subset in merged:
						merged[subset] += occur
					else:
						merged[subset] = occur
	else:
		with open(file, "rt") as f:
			for l in f:
				l = l.split()
				if len(l) == 4:
					subset = _subset(l[2], options)
					occur = float(l[3])
					if subset in merged:
						merged[subset] += occur
					else:
						merged[subset] = occur

	with open(outfile, "wt") as result:
		for subset in merged:
			result.write("{}\t{}\n".format(subset, merged[subset]))



if __name__ == "__main__":
	infile = "Pelagomonas_metaT-filters.tsv"
	suffix = "." + infile.split(".")[-1]
	outfile = infile.replace(suffix, "-merged" + suffix)
	mergeby(infile, outfile, "layer,station")