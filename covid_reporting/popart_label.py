import argparse
import os,re
import popart_label
from datetime import date, datetime
from Bio import SeqIO

# ideas: 	modify node sizes by including a sample count for each included sequence
# 			

def getIndex(list, num):
	#compares num to a list of ranges and returns the index of a match
    for idx, val in enumerate(list):
        if(num in range(val[0],val[1])):
            return idx
    return -1


def traitsheader(list):
	return """BEGIN TRAITS;
    Dimensions NTRAITS={};
    Format labels=yes missing=? separator=Comma; 
    TraitLabels {};
    Matrix
""".format(len(list), " ".join(list))

#note to self: use an "unknown" label for sequences without annotation

def networkheader(list):
	colorlist = ["#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c",
				 "#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5",
				 "#8c564b","#c49c94","#e377c2","#f7b6d2","#7f7f7f",
				 "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"]
	default_colors = {
	#lineages:
	"B.1": "#A6CEE3",
	"P.1": "#B2DF8A",
	"P.2": "#33A02C",
	"B.1.1": "#FFFF99",
	"B.1.1.7": "#F3444B",
	"B.1.1.153": "#B15928",
	"B.1.1.318": "#1F78B4",
	"B.1.160": "#FDBF6F",
	"B.1.177": "#FB9A99",
	"B.1.258": "#999999",
	"B.1.258.3": "#999999",
	"B.1.351": "#FF7F00",
	"B.1.525": "#C65911",
	"B.1.617": "#C65911",
	"B.1.617.2": "#CAB2D6",
	"AY.1": "#9053D1",
	"AY.2": "#9053D1",
	"AY.3": "#9053D1",
	"AY.4": "#128DA8",#
	"AY.43": "#128DA8",#
	"AY.5": "#1276A8",#
	"AY.6": "#9053D1",
	"AY.7.1": "#9053D1",
	"AY.8": "#9053D1",
	"AY.9": "#6A3D9A",
	"AY.10": "#9053D1",
	"AY.11": "#9053D1",
	"AY.12": "#1254A8",#
	"AY.20.1": "#D4FB79",
	"AY.21": "#1222A8",#
	"AY.CZ": "#6233F8",
	"AY.98": "#D4FB79",
	"AY.98.1": "#D4FB79",
	"AY.113": "#A6CEE3",
	"AY.121": "#B2DF8A",
	"AY.122": "#33A02C",
	"AY.125": "#D4FB79",
	"AY.126": "#FFFF99",
	"AY.127": "#F3444B",
	"AZ.2": "#338CC8",
	"BA.1": "#FF5DBC",
	"BA.1.1": "#FF5DBC",
	"BA.2": "#C90076",
	"C.36": "#ABA200",
	"C.36.3": "#DDD100",
	"XY": "#6233F8",
	#"Other": "#999999",
	#"None": "#BBBBBB",

	#regions:
	"JC": "#1f77b4",
	"JM": "#aec7e8", 
	"KA": "#ff7f0e",
	"KR": "#ffbb78",
	"LI": "#2ca02c",
	"MO": "#98df8a",
	"OL": "#d62728",
	"PA": "#ff9896",
	"PH": "#9467bd",
	"PL": "#c5b0d5",
	"ST": "#8c564b",
	"US": "#c49c94",
	"VY": "#e377c2",
	"ZL": "#f7b6d2",
	"neznamo": "#797979",
	"XY": "#FF00F0"
	}

	#get a region color or cycle through the base palette
	palette = [default_colors.get(x, colorlist[i%20]) for i,x in enumerate(list)]
	missing = ",".join([x for x in list if x not in default_colors])
	if missing != "":
		print("  missing colors for:", missing)
	#print(palette)
	string = """

Begin Network;
Dimensions ntax=666 nvertices=830 nedges=864 plotDim=0,0,4100,5600;
Format Font=Baskerville,10,-1,5,50,0,0,0,0,0 LegendFont=Baskerville,30,-1,5,50,0,0,0,0,0 VColour=#000000ff EColour=#000000ff BGColour=#00000000 VSize=30 EView=Dashes LPos=-1,1 LColours={};
""".format(",".join(palette)) 
#Font format: type, size, ?, ?, ?, italic, underlined, crossed, ?, ?
	
	return string


def reformat_date(date):
	yyyymmdd = r'(\d{4})-(\d{2})-(\d{2})'
	ddmmyyyy = r'(\d{2})/(\d{2})/(\d{4})'
	ddmmyyyyhyphen = r'(\d{2})-(\d{2})-(\d{4})'
	mmddyyyy = r'(\d{1,2})/(\d{1,2})/(\d{4})'
	ddmmyyyydotted = r'(\d{1,2})\.(\d{1,2})\.(\d{4})'

	try:
		if re.match(yyyymmdd, date):
			z = re.match(yyyymmdd, date)
			year = int(z.group(1))
			month = int(z.group(2))
			assert(month < 13), "Wrong format! {}".format(date)
			day = int(z.group(3))

		elif re.match(ddmmyyyy, date):
			z = re.match(ddmmyyyy, date)
			year = int(z.group(3))
			month = int(z.group(2))
			assert(month < 13), "Wrong format! {}".format(date)
			day = int(z.group(1))

		elif re.match(ddmmyyyyhyphen, date):
			z = re.match(ddmmyyyyhyphen, date)
			year = int(z.group(3))
			month = int(z.group(2))
			assert(month < 13), "Wrong format! {}".format(date)
			day = int(z.group(1))

		elif re.match(ddmmyyyydotted, date):
			z = re.match(ddmmyyyydotted, date)
			year = int(z.group(3))
			month = int(z.group(2))
			assert(month < 13), "Wrong format! {}".format(date)
			day = int(z.group(1))

		elif re.match(mmddyyyy, date):
			z = re.match(mmddyyyy, date)
			year = int(z.group(3))
			month = int(z.group(1))
			assert(month < 13), "Wrong format! {}".format(date)
			day = int(z.group(2))

		else:
			#print("Failed to parse date", date)
			return 0, 0, 0
	except:
		raise ValueError("\tdate format exception {}".format(date))

	consdate = year,month,day

	return consdate


def parse_pangolin(infile, infasta):
	seqset = {x.name for x in SeqIO.parse(infasta, "fasta")}
	lineage_d = {"neznamo": []}
	seen = set()

	with open(infile, "rt") as f:
		fset = [x for x in f.readlines() if len(x.split(",")) > 1]
		fdict = {x.split(",")[0]: x.split(",")[1].strip() for x in fset} #only last items will be stored
		for l in fdict:
			#l = l.strip().split(",")
			if l not in seqset:
				#print(l[0], "not in lineages.csv")
				continue
			seen.add(l)
			if fdict[l] == "None":
				fdict[l] = "neznamo"
			if fdict[l] not in lineage_d:
				lineage_d[fdict[l]] = []
			lineage_d[fdict[l]].append(l)
		unknown = list(seqset - seen)
		unknown.sort()
		lineage_d["neznamo"] += unknown

	print("PopARTLBL: sequences imported: {}".format(len(seqset))) 
	print("  Data found for {} sequences (out of {})".format(len(seen), len(seqset)))
	if len(seqset) - len(seen) > 0:
		print("  Missing: {}".format(",".join(seqset - seen)))
	print("PopARTLBL: lineages total {}".format(len(lineage_d)-1))
	return lineage_d


def parse_datafreeze(samplesdir, whichcol):
	"""Parses a set of XLSX in a folder to store information on sample collection
	date and site, potentially other data too.
	Columns:
	
	Accession ID		[0]
	GISAID_ID			[1]
	Collection date		[2]
	Location			[3]
	Gender				[4]
	Patient age			[5]
	Region code			[6]
	
	Args:
		samplesdir:		Directory with sample_sheet files.
			
	Returns:
		samples_dict:	Dictionary of accession->date.
	lineage_d = {}
	"""
	samples_dict = {}
	seen = 0
	_fileextension = ".tsv"
	_columns = ("id", "gisaid_id", "collection_date", "Location", "gender", "age", "rcode")

	files = [os.path.join(root, x) for root, dirs, files in os.walk(samplesdir)
			for x in files if _fileextension in x]
	
	for file in files:
		#check if file has exactly 20 columns:
		header = open(file, "rt").readline()
		if len(header.split("\t")) != 20:
			print("  ERROR! Incorrect number of columns! Skipping file: {}".format(file))
			continue
		else:
			print("Processing {}".format(file))
		with open(file, "rt") as f:
			samples = {l.split("\t")[0]: l.split("\t")[whichcol].strip() for l in f if len(l)>0}
			#trait-less sequences are improperly parsed by PopART!
			#this needs a different work-around:
			samples = {x: (y if y != "" else "neznamo") for x,y in samples.items()}
			samples = {x: y for x,y in samples.items() if y not in _columns}
			seen += len(samples.keys())
		
		if whichcol == 4: #this is a date column
			samples = {x: y for x,y in samples.items() if y != "neznamo"}
			for s, v in samples.items():
				#v is date in str
				y,m,d = reformat_date(v)
				try:
					date = datetime(y,m,d)
					s, y, w = str(s), date.strftime("%Y"), date.strftime("%V")
					if (y, w) not in samples_dict:
						samples_dict[(y, w)] = set()
					samples_dict[(y, w)].add(s)
				except ValueError:
					print("  ERROR! date format error, skipping: {}".format(v))
					continue
		else:
			for s, v in samples.items():
				s, v = str(s), str(v)
				if v not in samples_dict:
					samples_dict[v] = set()
				#s = s.replace("_", "-")
				#v = v.replace(",", ".")
				samples_dict[v].add(s)

	print("PopARTLBL: metadata found for {} sequences".format(seen))
	print("PopARTLBL: metadata items: {}".format(len(samples_dict)))
	return samples_dict


def parse_xls(samplesdir, whichcol):
	"""Parses a set of XLSX in a folder to store information on sample collection
	date and site, potentially other data too.
	valid = GISAID columns:
	Virus name / Isolate name
	Accession ID
	Type
	Clade
	Pango lineage
	Pango lineage version / date
	AA substitutions
	Variant
	Passage details/history
	ISIN identifier
	Receipt date
	Collection date
	Collector name
	Collecting institution
	Location
	Geographic latitude
	Geographic longitude
	ZIP code
	Sample capture status
	Host disease outcome
	Host
	Additional location information
	Sex
	Patient age
	Patient status
	Date of symptom onset
	Symptoms
	Clinical outcome
	Specimen source
	Travel history
	Cluster name
	Demographic patient info
	Clinical patient info
	Additional host information
	Sampling strategy
	Outbreak
	Last vaccinated
	Vaccine
	Treatment
	Motivation of sequencing
	CT value
	Sequencing technology
	Library preparation strategy
	Library preparation kit
	Assembly method, consensus generation
	Coverage
	Minimum sequencing depth to call sites
	FASTQ file name
	Duplicate sample
	Comment
	Originating lab
	Originating lab address
	Sample ID given by the originating laboratory
	Submitting lab
	Submitting lab address
	Sample ID given by the submitting laboratory
	Authors
	Submitter
	Submission Date
	Submitter address
	
	Args:
		samplesdir:		Directory with sample_sheet files.
			
	Returns:
		samples_dict:	Dictionary of accession->date.

	NOTE: only works with xlrd 1.2.0
	>> pip3 install xlrd==1.2.0 <<
	"""
	_filename = "sample_sheet" #could be also XLSX
	_sheetname = "Sheet1"
	_sample = "Sample ID given by the submitting laboratory"
	_location = "Location"
	_collection = "Collection date"
	usecols = [_sample, _location, _collection]
	if whichcol not in usecols:
		usecols.append(whichcol)

	files = [os.path.join(root, x) for root, dirs, files in os.walk(samplesdir)
			for x in files if _filename in x]
	files = [x for x in files if "~$" not in x] #drop temporary files

	samples_dict = {}
	for file in files:
		print(file)
		with open(file, "rb") as f:
			#https://en.wikipedia.org/wiki/List_of_file_signatures
			if f.read(2).decode("utf-8") != "PK":
				print(f.read(2))
				#raise ValueError("{} is not a XLSX file".format(file))
		try:
			samples = pd.read_excel(file,
					sheet_name=_sheetname,
					header=0,
					index_col=False,
					keep_default_na=True,
					usecols=usecols,
					)
		except IndexError:
			print("something went wrong, no idea what errors to expect")
		samples = samples[[_sample, whichcol]] #sample names must be first
		samples = dict(samples[samples[whichcol] != ""].values)
		
		if "date" in whichcol:
			for s, v in samples.items():
				#v is a Timestamp type, format to weeks
				s, y, w = str(s), v.strftime("%Y"), v.strftime("%V")
				samples_dict[s] = (y, w)
		else:
			for s, v in samples.items():
				s, v = str(s), str(v)
				#s = s.replace("_", "-")
				#v = v.replace(",", ".")
				samples_dict[s] = v
	
	return samples_dict


def clades(lineages_d, infasta, threshold):
	seqset = set(x.name for x in SeqIO.parse(infasta, "fasta"))

	lineages = [x for x in lineages_d.keys() if len(lineages_d[x]) >= threshold]
	lineages.sort()
	print("PopARTLBL: lineages filtered {} (>={} count)".format(len(lineages), threshold))
	if "neznamo" not in lineages:
		lineages.append("neznamo")

	# #remove empty datasets - not necessary for lineages which are derived from the fasta dataset
	# for l in lineages_d:
	# 	lineages_d[l] = [x for x in lineages_d[l] if x in seqset]

	# skipped = [x for x in lineages if len(lineages_d[x]) == 0]
	# lineages = [x for x in lineages if len(lineages_d[x]) > 0]
	# if skipped:
	# 	print("  No data found for {}".format(", ".join(skipped)))

	#another problem is with handling sequences that don't have a lineage assigned or fail the threshold
	#these need to be reported in the traits file as well
	other = [x for x in lineages_d.keys() if len(lineages_d[x]) < threshold]
	for lineage in other:
		if lineage == "neznamo":
			continue
		for seq in lineages_d[lineage]:
			lineages_d["neznamo"].append(seq)

	print("  sequences not passing threshold reassigned")
	#set headers for outputfiles
	#for nexus file - does not work but can be used to replace the traits block
	header_nex = traitsheader(lineages)
	tail_nex = ";\n\nEND;\n"

	#for csv file
	header_traits = ",{}\n".format(",".join(lineages))

	with open("popart_lineages.csv", "wt") as resultcsv,\
		open("popart_lineages.nex", "wt") as resultnex:
		resultcsv.write(header_traits)
		resultnex.write(header_nex)

		#print("writing taxa block by lineage")
		for lineage in lineages:
			for seq in lineages_d[lineage]:
				if seq not in seqset:
					continue
				traits = [lineages_d[x].count(seq) for x in lineages]
				traits = [str(x) for x in traits]
				resultcsv.write("{},{}\n".format(seq, ",".join(traits)))
				resultnex.write("{} {}\n".format(seq, ",".join(traits)))

		resultnex.write(tail_nex)
		resultnex.write(networkheader(lineages))


def dates(dates_d, infasta, freezetime):
	seqset = set(x.name for x in SeqIO.parse(infasta, "fasta"))
	seen = 0
	print("PopARTLBL: sequences imported: {}".format(len(seqset)))

	if freezetime == "today":
		thisy, thisw = int(date.today().strftime("%Y")), int(date.today().strftime("%V"))
		print("PopARTLBL: this week - {}".format(str(thisw)))
	else:
		y,m,d = reformat_date(freezetime)
		thisy, thisw = int(datetime(y,m,d).strftime("%Y")), int(datetime(y,m,d).strftime("%V"))
		print("PopARTLBL: this week {} based on setting {}".format(thisw, freezetime))

	if thisy == 2021:
		thisw += 53

	#ranges = ["1wa", "2wa", "3wa", "4wa", "2ma", "3ma", "4ma", "older"]
	ranges = [[thisw-1, thisw+1], 
			  [thisw-2, thisw-1], 
			  [thisw-3, thisw-2], 
			  [thisw-4, thisw-3],
			  [thisw-8, thisw-4], 
			  [thisw-12,thisw-8], 
			  [thisw-16,thisw-12], 
			  [0,thisw-15]]
	weeks = ["{}-{}".format((thisw-1) % 53,(thisw) % 53), 
			"{}".format((thisw-2) % 53), 
			"{}".format((thisw-3) % 53), 
			"{}".format((thisw-4) % 53), 
			"{}-{}".format((thisw-8) % 53,(thisw-5) % 53), 
			"{}-{}".format((thisw-12) % 53,(thisw-9) % 53), 
			"{}-{}".format((thisw-16) % 53,(thisw-11) % 53), 
			"older"]

	adjusted_d = {x: [] for x in weeks}

	#populate adjusted weeks:
	for (y,w) in dates_d:
		items = dates_d[(y,w)]
		seen += len([x for x in items if x in seqset])
		year, week = int(y), int(w)
		if year == 2021 and week != 53:
			week += 53
		if week > thisw:
			weirditems = [x for x in items if x in seqset]
			if len(weirditems) > 0:
				print("  Suspicious date Y{}-W{}, items: {}".format(y, w, weirditems))
			else:
				print("  Suspicious date Y{}-W{} in metadata, however all items missing from input fasta".format(y, w))
		weekstr = weeks[getIndex(ranges, week)]
		adjusted_d[weekstr].extend(items)

	#remove empty datasets
	for w in adjusted_d:
		adjusted_d[w] = [x for x in adjusted_d[w] if x in seqset]

	skipped = [x for x in weeks if len(adjusted_d[x]) == 0]
	weeks = [x for x in weeks if len(adjusted_d[x]) > 0]
	if skipped:
		print("  No seqs found for {}".format(", ".join(skipped)))

	print("  Data found for {} sequences (out of {})".format(seen, len(seqset)))
	if len(seqset) - len(seen) > 0:
		print("  Missing: {}".format(",".join(seqset - seen)))
	#set headers for outputfiles
	#for nexus file - does not work but can be used to replace the traits block
	header_nex = traitsheader(weeks)
	tail_nex = ";\n\nEND;\n"

	#for csv
	header_traits = ",{}\n".format(",".join(weeks))

	with open("popart_dates.csv", "wt") as resultcsv,\
		open("popart_dates.nex", "wt") as resultnex:
		resultcsv.write(header_traits)
		resultnex.write(header_nex)

		for weekstr in weeks:
			for seq in adjusted_d[weekstr]:
				if seq not in seqset:
					continue
				#print(seq)
				traits = [adjusted_d[x].count(seq) for x in weeks]
				traits = [str(x) for x in traits]
				resultcsv.write("{},{}\n".format(seq, ",".join(traits)))
				resultnex.write("{} {}\n".format(seq, ",".join(traits)))

		resultnex.write(tail_nex)


def regions(regions_d, infasta):
	seqset = set(x.name for x in SeqIO.parse(infasta, "fasta"))
	seen = set()
	print("PopARTLBL: sequences imported: {}".format(len(seqset)))

	regions = [x for x in regions_d.keys()]
	regions.sort()

	#remove empty datasets
	for r in regions_d:
		regions_d[r] = [x for x in regions_d[r] if x in seqset]

	skipped = [x for x in regions if len(regions_d[x]) == 0]
	regions = [x for x in regions if len(regions_d[x]) > 0]
	if skipped:
		print("  No seqs found for {}".format(", ".join(skipped)))

	#set headers for outputfiles
	#for nexus file - does not work but can be used to replace the traits block
	header_nex = traitsheader(regions)
	tail_nex = ";\n\nEND;\n"

	#for csv file
	header_traits = ",{}\n".format(",".join(regions))

	with open("popart_regions.csv", "wt") as resultcsv,\
		open("popart_regions.nex", "wt") as resultnex:
		resultcsv.write(header_traits)
		resultnex.write(header_nex)

		for region in regions:
			for seq in regions_d[region]:
				if seq not in seqset:
					continue
				#print(seq)
				if seq in seen:
					continue
				seen.add(seq)
				traits = [regions_d[x].count(seq) for x in regions]
				traits = [str(x) for x in traits]
				resultcsv.write("{},{}\n".format(seq, ",".join(traits)))
				resultnex.write("{} {}\n".format(seq, ",".join(traits)))

		resultnex.write(tail_nex)
		resultnex.write(networkheader(regions))
	print("  Data found for {} sequences (out of {})".format(len(seen), len(seqset)))
	if len(seqset) - len(seen) > 0:
		print("  Missing: {}".format(",".join(seqset - seen)))


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-l', '--lineage', help='Lineage as PopART trait', action='store_true')
	group.add_argument('-r', '--region', help='Region as PopART trait', action='store_true')
	group.add_argument('-d', '--date', help='Date/week as PopART trait', action='store_true')

	parser.add_argument('-f', '--fastafile', help='Fasta as input (for filtering)', required=True)
	parser.add_argument('--lineagefile', help='Lineage report as input', default="lineage_report.csv")
	parser.add_argument('--freezedate', help='Lineage report as input', default="today")
	parser.add_argument('--metadatadir', help='Sample sheet(s) directory', default=".")
	#parser.add_argument('-i', '--infile', help='Nexus infile', default="")
	parser.add_argument('-t', '--threshold', help='Threshold for lineage inclusion', default=2)

	args = parser.parse_args()
	threshold = int(args.threshold)

	global ignore
	ignore = [""]

	if args.lineage:
		print("PopARTLBL: lineage traits to parse")
		clades(parse_pangolin(args.lineagefile, args.fastafile), args.fastafile, threshold)
	elif args.region:
		print("PopARTLBL: region traits to parse")
		regions(parse_datafreeze(args.metadatadir, 10), args.fastafile)
	# 	regions(parse_xls(args.metadatadir, "Location"))
	elif args.date:
		print("PopARTLBL: date traits to parse relative to {}".format(args.freezedate))
		dates(parse_datafreeze(args.metadatadir, 4), args.fastafile, args.freezedate)
	#	dates(parse_xls(args.metadatadir, "Collection date"))


if __name__ == '__main__':
	main()

