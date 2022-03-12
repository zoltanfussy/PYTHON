import os, argparse, sys, re
import pandas as pd
from Bio import SeqIO,AlignIO
from datetime import date, datetime, timedelta
import covid_phylo
#from pprint import pprint as pp
#xlrd required! pip3 install xlrd==1.2.0


def print_and_log(string):
	with open("_run.log", "at") as result:
		result.write("{}\n".format(string))
	print(string) 

def _reformat_date(date):
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
			#print("  Failed to parse date {}".format(date))
			return 2020, 1, 1 #before the sequencing started
	except:
		raise ValueError("\tdate format exception {}".format(date))

	consdate = year,month,day

	return consdate


def _parse_dates(dates_csv, col_id, separator=",", substring_date=True):
	"""Parses a CSV having AC in the first and collection_date
	in the last column given by col_id.
	
	Args:
		dates_csv:	CSV filename.
		col_id: 		Column ID of the date
			
	Returns:
			Dictionary of accession->date.
	"""
	dates_dict = {}
	with open(dates_csv) as f:
		for l in f:
			#Accession,Species,Geo_Location,Host,Isolation_Source,Collection_Date
			l = l.strip().split(separator)
			if len(l) < col_id + 1:
				continue
			
			seqid = l[0]
			if substring_date:
				# Take just a part of the date string matching YY-MM-DD 
				# => GISAID fasta header format
				coldate = l[col_id][-8:]
			else:
				coldate = l[col_id]
			dates_dict[seqid] = coldate
	
	return dates_dict
	
		
def _parse_csv(dates_csv, col_id, separator=","):
	"""Parses a CSV having AC in the first and collection_date
	in the last column given by col_id.
	
	Args:
		dates_csv:		CSV filename.
		col_id: 		Column ID of the date
			
	Returns:
			Dictionary of accession->date.
	"""
	dates_dict = {}
	with open(dates_csv) as f:
		for l in f:
			#Accession,Species,Geo_Location,Host,Isolation_Source,Collection_Date
			l = l.strip().split(separator)
			if len(l) < col_id + 1:
				continue
			
			seqid = l[0]
			coldate = l[col_id]
			dates_dict[seqid] = coldate
	
	return dates_dict
	
		
def _pandas_csv(dates_csv, col_id="collection_date", separator=","):
	"""Parses a CSV having AC in the first and collection_date
	in the last column given by col_id.
	
	Args:
		dates_csv:		CSV filename.
		col_id: 		Column ID of the date
			
	Returns:
			Dictionary of accession->date.
	"""
	df = pd.read_csv(dates_csv, delimiter=separator)
	df = df[["fasta_id",col_id]]
	df.set_index("fasta_id")
	dates_dict = dict(df.values)
	
	return dates_dict
	
		
def _parse_xls(samplesdir, whichcol):
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
		print_and_log(file)
		with open(file, "rb") as f:
			#https://en.wikipedia.org/wiki/List_of_file_signatures
			if f.read(2).decode("utf-8") != "PK":
				print_and_log(f.read(2))
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
			print_and_log("something went wrong, no idea what errors to expect")
		samples = samples[[_sample, whichcol]] #sample names must be first
		samples = dict(samples[samples[whichcol] != ""].values)
		
		if "date" in whichcol:
			for s, v in samples.items():
				#v is a Timestamp type, format to YY-MM-DD
				s, v = str(s), v.strftime("%y-%m-%d")
				samples_dict[s] = v
		else:
			for s, v in samples.items():
				s, v = str(s), str(v)
				#s = s.replace("_", "-")
				#v = v.replace(",", ".")
				samples_dict[s] = v
	
	return samples_dict

		
def _parse_cogcz_seqname(seqname):
	"""Modifies assembly sequence name to a shorter one,
	retaining information on AC, region of origin, and collection date.
	
	Args:
			Seqname.
			
	Returns:
			Modified seqname.
	
	"""
	#datesdb and regiondb inherited from cz_regions
	try:
		seqid = seqname
		coldate = datesdb.get(seqid, "")
		region = regiondb.get(seqid, "CZE")
		if coldate == "":
			seqname = "{}@{}".format(seqid,region)
		else:
			seqname = "{}@{}@{}".format(seqid,region,coldate)
	except:
		#print("could not parse", seqname)
		print_and_log("WARNING: Could not parse {}".format(seqname))

	return seqname
	
	
def _parse_gisaid_seqname(seqname):
	"""Modifies sequence name to a shorter one,
	retaining information on AC, country of origin, and collection date.
	
	Args:
			Seqname.
			
	Returns:
			Modified seqname.
	
	"""
	seqname = seqname.replace(" ", "_")
	try:
		country = seqname.split("/")[1]
		seqid = seqname.split("|")[-2]
		coldate = seqname[-8:]
		return "{}@{}@{}".format(seqid,country,coldate)
	except IndexError:
		#print("could not parse", seqname)
		seqidv = seqname.split("_")[0]
		seqid = seqidv.split(".")[0]
		if seqid in datesdb:
			coldate = datesdb[seqid]
			country = "_".join(seqname.split("_")[1:])
			seqname = "{}@{}@{}".format(seqid,country,coldate)
		elif "NC_045512" in seqname:
			pass
		else:
			print_and_log("WARNING: Could not parse {}".format(seqname))
			seqname = seqname.replace(seqidv + "_", seqidv + "@")
		return seqname
	
	
def filter(infile, outfile, filterfile):
	"""Deduplicates an infile Fasta, using a filter
	list to write unwanted results to a different file.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
		filterfile:	File containing sequences to skip
		
	Returns:
		N/A 		Writes two files
	
	"""
	seqset = set()
	seen = set()
	with open (filterfile) as f:
		filtset = {x.strip().split("\t")[0] for x in f.readlines()}
	#print(filtset)
	print_and_log("PHYLO: Filter infasta: {}, remove {} sequences.".format(infile, len(filtset)))
	with open(outfile, "w") as result,\
	open("f" + outfile, "w") as result2:
		c, f = 0, 0
		for seq in SeqIO.parse(infile, "fasta"):
			seqname = seq.name
			if  seqname  in seqset:
				c += 1
				continue
			if  seqname  in filtset:
				f += 1
				seen.add(seqname)
				result2.write(">{}\n{}\n".format(seqname, seq.seq))
				continue
			if len(seqname) >50:
				print_and_log("Seq too long for BLAST! {}".format(seqname))
			result.write(">{}\n{}\n".format(seqname, seq.seq))
			seqset.add(seqname)
	print_and_log("  {} sequences filtered, {} sequences dereplicated".format(f, c))
	if len(filtset - seen) > 0:
		print_and_log("Some sequence names not seen: {}".format(filtset - seen))


def trim(infile, outfile, seqn_threshold, datafile, clean_gaps, verbose):
	"""Trims an alignment Fasta, using a filter
	list to trash unwanted results. Also, sequences
	with high number of unknowns are removed.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes two files
	
	"""
	print_and_log("PHYLO: Trimming\nProcess alignment: {}".format(infile))
	#for testing only, this writes a summary of failed sequences to a reportfile
	#datafile = "df-20210506-v2.tsv"
	if os.path.exists(datafile):
		#datesdb = _parse_csv(datafile, col_id=2, separator=separator) #column order before 21.5.
		datesdb = _parse_csv(datafile, col_id=4, separator=separator) #after 21.5.!
		#regiondb = _parse_csv(datafile, col_id=-1, separator=separator) #rcode always last column
		regiondb = _parse_csv(datafile, col_id=-2, separator=separator) #after 21.5.!
	else:
		datesdb = {}
		regiondb = {}
	
	gapchar = "-" if clean_gaps else ""
	seqset = set()
	#print(filtset)
	alignment = AlignIO.read(infile, "fasta")
	if reflength == 29823:
		# The reference genome has 29903 nt.
		# By default, this function uses a nextalign input and the ends need to be trimmed
		# because of poor quality but this can be disabled here for manual alignments 
		# by passing a different reflength using --trim_threshold
		trim = alignment[:, 40:-40]
	else:
		trim = alignment[:, :]
	with open(outfile, "w") as result,\
		open("trim_result.log", "at") as report:
		report.write("Trimming started {}\n".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
		p,d,c = 0, 0, 0
		seqlen_threshold = int(por_n*reflength)
		#seqn_threshold = 30 # for relative threshold: round((1-por_n)*reflength)
		print_and_log("  Min len {}, max Ns {}".format(seqlen_threshold, seqn_threshold))
		for seq in trim:
			seqname = seq.description.replace(" ", "_")
			if seqname in seqset: #duplicate
				c += 1
				continue
			if len(seq.seq) < seqlen_threshold:
				#genome length not met
				report.write("{}\t{}\t{}\t{}: {}\n".format(seqname, 
														   datesdb.get(seqname, ""), 
														   regiondb.get(seqname, ""), 
														   "sequence too short",
														   len(seq.seq)
														   ))
				d += 1
				continue
			#if seq.seq.count("N") > seqn_threshold:
			if sum(seq.seq.upper().count(x) for x in "NWSMKRY-") > seqn_threshold:
				#too many N's
				#print("{}: {} knowns, {} unknowns".format(len(seq.seq), 
				#										   sum(seq.seq.count(x) for x in "ACTG-"), 
				#										   sum(seq.seq.count(x) for x in "NWSMKRY")))
				report.write("{}\t{}\t{}\t{}: {}\n".format(seqname, 
														   datesdb.get(seqname, ""), 
														   regiondb.get(seqname, ""), 
														   "too many unknowns",
														   sum(seq.seq.upper().count(x) for x in "NWSMKRY")
														   ))
				d += 1
				continue

			#print("{}: {} knowns, {} unknowns".format(len(seq.seq), 
			#										  sum(seq.seq.upper().count(x) for x in "ACTG-"), 
			#										  sum(seq.seq.upper().count(x) for x in "NWSMKRY")))
			result.write(">{}\n{}\n".format(seqname, str(seq.seq).replace(gapchar, "")))
			report.write("{}\t{}\t{}\t{}\n".format(seqname, 
												   datesdb.get(seqname, ""), 
												   regiondb.get(seqname, ""),
												   "sequence passed"
												   ))
			seqset.add(seqname)
			p += 1
		print_and_log("  {0} sequences failed QC at {1}% seq length and {2} Ns ({3} passed)".format(d, por_n, seqn_threshold, p))
	#AlignIO.write(trim, outfile, "fasta")
	if verbose:
		for seq in trim:
			if len(seq.seq) < seqlen_threshold:
				#sequence too short
				print("{} too short".format(seq.name))
			elif sum(seq.seq.upper().count(x) for x in "NWSMKRY") > seqn_threshold:
				#print("{}: {} knowns, {} unknowns".format(len(seq.seq), sum(seq.seq.count(x) for x in "ACTG-"), sum(seq.seq.count(x) for x in "NWSMKRY")))
				#too many N's
				print("{} too many ambigs".format(seq.name))
				continue
			else:
				print("{} OK".format(seq.name))	


def sort_by_date(infile, freezedate, span, datafile):
	"""Deduplicates an infile Fasta, using a filter
	list to trash unwanted results. Also, sequences
	with high number of unknowns are removed.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Sort recent sequences on collection date.")
	#for testing only:
	#datafile = "df-20210506-v2.tsv"
	if os.path.exists(datafile):
		datesdb = _parse_csv(datafile, col_id=4, separator=separator)
		datesdb = {x: _reformat_date(y) for x,y in datesdb.items()} # Reformat datesdb to unified format
		datesdb = {x: datetime(y[0], y[1], y[2]).date() for x,y in datesdb.items()} # Convert to datetime object
	else:
		quit("Metadata file missing")
	
	if freezedate == "":
		freezedate = date.today()
	else:
		y,m,d = _reformat_date(freezedate)
		freezedate = datetime(y,m,d).date()
	lowerdate = freezedate - timedelta(days=span)
	nodate = datetime(2020,1,1).date()


	seqset = set()
	#print(filtset)
	print_and_log("PHYLO: Sorting by date\nProcessing file: {}".format(infile))
	with open("czcurrent.tmp.fasta", "wt") as latest,\
		 open("czformer.tmp.fasta", "wt") as earlier,\
		 open("sort_result.log", "at") as report:
		p,d,n,e = 0, 0, 0, 0 #passed, duplicates, no_data, earlier
		print_and_log("  Timeframe {} - {}".format(lowerdate, freezedate))
		for seq in SeqIO.parse(infile, "fasta"):
			#assuming seqnames are matching the metadata!
			#seqname = seq.description.replace(" ", "_")
			seqname = seq.name
			seqdate = datesdb.get(seqname, nodate)
			if seqname in seqset: #duplicate
				d += 1
				continue
			if seqdate > freezedate:
				#sequence from future?
				earlier.write(">{}\n{}\n".format(seqname, str(seq.seq)))
				report.write("{}\t{}: {}\n".format(seqname, 
												   seqdate, 
												   "sequence date misformatted"
												   ))
				n += 1
			elif seqdate == nodate:
				earlier.write(">{}\n{}\n".format(seqname, str(seq.seq)))
				report.write("{}\t{}: {}\n".format(seqname, 
												   seqdate, 
												   "sequence date misformatted"
												   ))
				n += 1
			elif seqdate >= lowerdate:
				#sequence from the timespan
				latest.write(">{}\n{}\n".format(seqname, str(seq.seq)))
				report.write("{}\t{}: {}\n".format(seqname, 
												   seqdate, 
												   "sequence within set timeframe"
												   ))
				p += 1
			else:
				#sequence from earlier, or date unknown
				earlier.write(">{}\n{}\n".format(seqname, str(seq.seq)))
				report.write("{}\t{}: {}\n".format(seqname, 
												   seqdate, 
												   "sequence from earlier time"
												   ))
				e += 1
			seqset.add(seqname)

		print_and_log("  {0} sequences failed date check, {1} collected earlier, {2} passed.".format(n, e, p))


def representatives(infile, outfile, repsfile):
	"""Processes an infile Fasta - renames 
	sequence names to show a number of sequences 
	assigned to the cluster representative.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Process representatives list")
	
	with open(repsfile) as f:
		repslist = [x.split("\t")[0] for x in f.readlines()]
		repsset = set(repslist)
		repsdb = {x: repslist.count(x) for x in repsset}
	#pp(repsdb)
	print_and_log("PHYLO: Found {} representatives".format(len(repsdb)))

	with open(outfile, "w") as result:
		for seq in SeqIO.parse(infile, "fasta"):
			count = repsdb[seq.name]
			seqname = "{}_x{}".format(seq.name, count)
			if len(seqname) >50:
				print_and_log("WARNING: Seq too long for BLAST! {}".format(seqname))
			result.write(">{}\n{}\n".format(seqname, seq.seq))


def world_regions(infile, outfile, seqn_threshold, datafile, col_id):
	"""Processes an infile Fasta - deduplicates and renames 
	sequence names to have a unified ID@REGION@DATE format.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Extract date information for GenBank accessions")
	if col_id == "":
	#column id different than -1!
		col_id = -1

	if datafile == "":
	#column id different than -1!
		datafile = r"GenBank_data.csv" #obsolete

	global datesdb
	datesdb = {}#_parse_dates(datafile, col_id)
	#print_and_log("PHYLO: Imported {} items".format(len(datesdb)))
	
	print_and_log("PHYLO: Process infasta: {}".format(infile))
	seqset = set()
	with open(outfile, "w") as result:
		p,d = 0, 0
		seqlen_threshold = int(por_n*reflength)
		#seqn_threshold = 30 #round((1-por_n)*reflength)
		print_and_log("  Min len {}, max Ns {}".format(seqlen_threshold, seqn_threshold))
		for seq in SeqIO.parse(infile, "fasta"):
			if len(seq.seq) < seqlen_threshold:
				#sequence too short
				d += 1
				continue
			#if seq.seq.count("N") > seqn_threshold:
			if sum(seq.seq.upper().count(x) for x in "NWSMKRY") > seqn_threshold:
				#print("{}: {} knowns, {} unknowns".format(len(seq.seq), sum(seq.seq.upper().count(x) for x in "ACTG-"), sum(seq.seq.upper().count(x) for x in "NWSMKRY")))
				#too many N's
				d += 1
				continue
			if "@" not in seq.name:
				seqname = _parse_gisaid_seqname(seq.description)
			else:
				seqname = seq.name
			if len(seqname) >50:
				print_and_log("WARNING: Seq too long for BLAST! {}".format(seqname))
			if not seqname in seqset:
				result.write(">{}\n{}\n".format(seqname, seq.seq))
			p += 1
			seqset.add(seqname)
		print_and_log("  {0} sequences failed QC at {1}% seq length and {2} Ns ({3} passed)".format(d, por_n, seqn_threshold, p))


def cz_regions(infile, outfile, seqn_threshold, datafile, verbose):
	"""Processes an infile Fasta - renames 
	sequence names to have a unified ID@REGION@DATE format.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Extract collection date and location for local accessions")

	if datafile == "":
		#this was meant to find all xlsx files in a folder and its subfolders
		#but this approach was discontinued
		#set proper file path!
		datadir = r"/storage/vestec1-elixir/projects/cogcz/seq"
		datadir = "."

	global datesdb
	#datesdb = _parse_xls(datadir, "Receipt date")
	datesdb = _parse_csv(datafile, col_id=4, separator=separator)
	datesdb = {x: y.replace("/","-") for x,y in datesdb.items()}
	print_and_log("PHYLO: Imported {} items".format(len(datesdb)))

	global regiondb
	#regiondb = _parse_xls(datadir, "Location")
	regiondb = _parse_csv(datafile, col_id=-2, separator=separator)
	regiondb = {(x if True else x): ("CZE-"+y if y != "" else "CZE") for x,y in regiondb.items()}

	print_and_log("PHYLO: Process infasta: {}".format(infile))
	with open(outfile, "w") as result:
		p,d = 0, 0
		seqlen_threshold = int(por_n*reflength)
		#seqn_threshold = 30 #round((1-por_n)*reflength)
		print_and_log("  Min len {}, max Ns {}".format(seqlen_threshold, seqn_threshold))
		for seq in SeqIO.parse(infile, "fasta"):
			if len(seq.seq) < seqlen_threshold:
				#sequence too short
				d += 1
				continue
			#if seq.seq.count("N") > seqn_threshold:
			if sum(seq.seq.count(x) for x in "NWSMKRY") > seqn_threshold:
				#print("{}: {} knowns, {} unknowns".format(len(seq.seq), sum(seq.seq.count(x) for x in "ACTG-"), sum(seq.seq.count(x) for x in "NWSMKRY")))
				#too many N's
				d += 1
				continue
			if "@" not in seq.name:
				seqname = _parse_cogcz_seqname(seq.description)
			else:
				seqname = seq.name
			if len(seqname) >50:
				print_and_log("WARNING: Seq too long for BLAST! {}".format(seqname))
			p += 1
			result.write(">{}\n{}\n".format(seqname, seq.seq))
		print_and_log("  {0} sequences failed QC at {1}% seq length and {2} Ns ({3} passed)".format(d, por_n, seqn_threshold, p))
	if verbose:
		for seq in SeqIO.parse(infile, "fasta"):
			if len(seq.seq) < seqlen_threshold:
				#sequence too short
				print("{} too short".format(seq.name))
			elif sum(seq.seq.count(x) for x in "NWSMKRY") > seqn_threshold:
				#print("{}: {} knowns, {} unknowns".format(len(seq.seq), sum(seq.seq.count(x) for x in "ACTG-"), sum(seq.seq.count(x) for x in "NWSMKRY")))
				#too many N's
				print("{} too many ambigs".format(seq.name))
				continue
			else:
				print("{} OK".format(seq.name))	
	

def mark_treelabel(infile, outfile, modsfile, col_id=1, separator=","):
	"""Processes an infile Fasta - renames 
	sequence names to show a number of sequences 
	assigned to the cluster representative.
	
	Args:
		infile:		Input tree name
		outfile: 	Output tree name
		modsfile:	Modifications as table
		separator:	Modifications table separator
		col_id:		Index of the column to be used as modifier
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Process label modification file")
	global modsdb

	modsdb = _parse_csv(modsfile, col_id=col_id, separator=separator)
	modsdb = {x.replace("@", "_"): y for x,y in modsdb.items()}

	#pp(modsdb)
	print_and_log("PHYLO: Found {} representatives".format(len(modsdb)))
	with open(infile) as f, open(outfile, "w") as result:
		outtree = f.read()
		for label, modifier in modsdb.items():
			newlabel = "{}_{}".format(label, modifier)
			outtree = outtree.replace(label, newlabel)
		result.write(outtree)


def dereplicate(infile, outfile):
	"""Deduplicates an infile Fasta on seqnames.
	
	Args:
		infile:		Input fasta name
		outfile: 	Output fasta name
			
	Returns:
		N/A 		Writes a file
	
	"""
	print_and_log("PHYLO: Derep infasta: {}".format(infile))
	seqset = set()
	with open(outfile, "w") as result:
		c = 0
		for seq in SeqIO.parse(infile, "fasta"):
			seqname = seq.description.replace(" ", "_")
			#if len(seqname) >50:
			#	print("Seq too long for BLAST! {}".format(seqname))
			if seqname in seqset:
				c += 1
				continue
			result.write(">{}\n{}\n".format(seqname, seq.seq))
			seqset.add(seqname)
	print_and_log("  {} sequences dereplicated".format(c))


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-i', '--infasta', help='Fasta to be processed', required=True)
	parser.add_argument('-o', '--outfasta', help='Output filename', default='')
	parser.add_argument('-d', '--data', help='Data file', default='')
	parser.add_argument('--today', help='Datafreeze date', default='')
	parser.add_argument('--timespan', help='Timeframe in days', default=30)
	parser.add_argument('--separator', help='Separator', default='\t')
	parser.add_argument('--column', help='Dates column in data file', default='')
	parser.add_argument('--trim_threshold', help='trimming threshold for Ns in sequences', default=30)
	parser.add_argument('--ali_length', help='alignment_length', default=29823)
	parser.add_argument('--clean_gaps', help='Clean gaps after trimming', action='store_true')
	parser.add_argument('-v', '--verbose', help='Verbose mode', action='store_true')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-f', '--filter', help='Filtering list', default='')
	group.add_argument('-t', '--trim', help='Perform trim', action='store_true')
	group.add_argument('-s', '--sort_by_date', help='Sort input fasta by date', action='store_true')
	group.add_argument('-r', '--representatives', help='Representatives list by mmseqs', default='')
	group.add_argument('-w', '--world_regions', help='Process data with world regions', action='store_true')
	group.add_argument('-cz', '--cz_regions', help='Process data with czech regions', action='store_true')
	group.add_argument('-m', '--mark_treelabel', help='Modify tree labels', action='store_true')

	args = parser.parse_args()

	infile = args.infasta

	if args.outfasta == "":
		outfile = args.infasta + ".uniq"
	else:
		outfile = args.outfasta

	if args.column != "":
		col_id = int(args.column)
	else:
		col_id = "" #reverting to each function's default

	datafile = args.data

	global por_n
	por_n = 0.995

	global reflength
	reflength = int(args.ali_length) #after trimming the first and last 40 bases!

	global separator
	separator = args.separator

	#print("Note that -r, -m, -cz and -w are mutually exclusive!")

	if args.representatives != "":
		#implicitly adds representatives
		representatives(infile=infile, outfile=outfile, repsfile=args.representatives)
	elif args.filter != "":
		#implicitly adds representatives
		filter(infile=infile, outfile=outfile, filterfile=args.filter)
	elif args.sort_by_date:
		sort_by_date(infile=infile, freezedate=args.today, span=int(args.timespan), datafile=datafile)
	elif args.trim:
		trim(infile=infile, outfile=outfile, seqn_threshold=int(args.trim_threshold), datafile=datafile, clean_gaps=args.clean_gaps, verbose=args.verbose)
	elif args.cz_regions:
		cz_regions(infile=infile, outfile=outfile, seqn_threshold=int(args.trim_threshold), datafile=datafile, verbose=args.verbose)
	elif args.world_regions:
		world_regions(infile=infile, outfile=outfile, seqn_threshold=int(args.trim_threshold), datafile=datafile, col_id=col_id)
	elif args.mark_treelabel:
		mark_treelabel(infile=infile, outfile=outfile, modsfile=datafile, separator=separator) #col_id=-1
	else:
		print_and_log("None of  -r / -cz / -w / -m chosen, dereplicating")
		dereplicate(infile=infile, outfile=outfile)

if __name__ == "__main__":
	main()

