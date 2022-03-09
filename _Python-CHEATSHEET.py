type(some_object)
if isinstance(function, str): #determine if instance is string/list/set/whatever
#The difference between type() and isinstance() is that isinstance() returns True even for instances of subclasses that inherit the class specified in the second argument.
#i.e., isinstance will be true for both isinstance(True, bool) and isinstance(True, int), since boolean values are subclass of int
print(os.path.isdir("/home/el"))
print(os.path.exists("/home/el/myfile.txt"))
path, file = os.path.split(filepath)

#non-ascii character error:
import io
notincluded = io.open("abstracts.txt", "r", encoding="utf-8") #this is an analog of open(file).read()

import csv
csv.reader(open("cafe_table.tsv"), delimiter="\t", skipinitialspace=True)
#or read csv/tsv as dictionary:
with open('enz2annot.txt', 'r') as f:
	reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
	next(reader) #another way to skip first line
	enz2annot = {r[0]: r[1] for r in reader}

#import matplotlib without Xcode library
import matplotlib
matplotlib.use("Agg")

#read csv with numpy:
import numpy as np
path = 'data/population.csv'
data = np.genfromtxt(path, delimiter=',', names=True)

#read csv with pandas:
import pandas as pd
path = 'data/population.csv'
df = pd.read_csv(path, skiprows=1, delimiter="\t", names=["freq", "count"])
# Set the country code as index of the DataFrame
df = df.set_index('Country Code')
df = df.loc[:,['Jan','Feb', 'Mar']] # use to rearrange the layer ordering

#>> melt is to merge several value columns into one, yielding two categories and one value column
melted = pd.melt(df, "layer", var_name="value") 

a_list.append() vs. a_list.extend(["another", "list"]) #append adds whatever argument as a single item

#do stuff for pandas subsets!
for i,k in enumerate(df.column.unique())
    groups = df.groupby('column')
    for name, group in groups:
        print(name, group) #group is a df subset!

break vs continue #break ends the for loop, continue ends the cycle for the current item
pass #use with if statements

import re
IPSpattern = r'"InterPro: (\w+)'
taxonpattern = r'\[(.+)\]' #to find taxons in fasta descriptions
genbankpattern = r'[A-Z]{1,3}(_)?\d+(\.\d)?'
taxon = re.search(taxonpattern, desc).group(1)
if taxon:
	taxon = "_".join(taxon.split()[:2])
	desc = f'{taxon}_{desc.split(" [")[0]}'
#this is to find the matching string
IPSID = re.search(IPSpattern, somestring).group(1)
#this is to iterate though matches
for hit in re.findall(IPSpattern, somestring):
	blabla
#POZOR, re.match hleda na zacatku stringu

remove() vs discard()

#zipping two lists into a list of tuples:
gradebook = list(zip(libraries, completion))

import os
#to run os commands
if os.path.isfile("{}.nhr".format(args.database)) == False:
	os.system('makeblastdb -in {} -out {} -dbtype nucl -parse_seqids'.format(args.database, args.database))
else:
	print "skipping makeblastdb... database exists"

cmd = 'blastn'
db = 'VBRAnt'
query = 'CryptoDB-34_VbrassicaformisCCMP3155_AnnotatedProteins.fasta'
out = 'vbra_all_blast_mmetsp.xml'
evalue = 10
outfmt = 5
word_size = 4
threads = 4
maxevalue = 0.01
maxhits = ''
if os.path.isfile("vbra_all_blast_mmetsp.xml") == False:
	print('{} -query {} -db {} -out {} -evalue {} -outfmt {} -word_size {} -num_threads {}'.format(cmd, query, db, out, evalue, outfmt, word_size, threads))
	#or
	os.system(f'{cmd} -query {query} -db {db} -out {out} -evalue {evalue} -outfmt {outfmt} -word_size {word_size} -num_threads {threads}')

#or
os.system('targetp -P {0}-out.fasta > {0}-targetp.txt'.format(prefix))

os.chdir("some/dir")
files = os.listdir('.')
files = [f for f in files if f.endswith(".fasta")]

#home = "/Users/zoliq/ownCloud/"
home = "/Volumes/zoliq data/OwnCloud/"
wd = home + "genomes/chromera/plastid proteome/"
os.chdir(wd)

if os.stat(fname).st_size > 0: #check for file size
	something()

if any(x in rank for x in goodgroups)

directory_walker = os.walk('.')
for d in directory:
	#d[0] current directory string, d[1] subdirectories list, d[2] current directory file list
	print(d)

x = "NNNcagat"
prve_N = x.find('N')
posledne_N = x.rfind('N')
print(prve_N, posledne_N) #indexy


import argparse
#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-o', '--outfile', help='Output filename', default="")
parser.add_argument('-b', '--bootstrap', help='Boostrap calculation', action='store_true') # to only look for presence/absence of this argument
group = parser.add_mutually_exclusive_group()
group.add_argument('-f', '--filter', help='Filtering list', default='')
group.add_argument('-t', '--trim', help='Perform trim', action='store_true')


args = parser.parse_args()

infile = args.infile
if args.outfile == "":
	suffix = infile.split(".")[-1]
	outfile = infile.replace(suffix, "renamed.fasta")
num_seqs = args.num_seqs

##################################
#### Create working directory ####
##################################

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
if args.directory == ".":
	print("changing to default directory")
	defdir = "Robolab/phylo/"
	wd = home + defdir
	os.chdir(wd)
else:
	os.chdir(args.directory)


import sys
#communicating with the system
#another way of importing arguments from command line
with open(sys.argv[1]) as infile, open(sys.argv[2]) as outfile:
	something
#determine if operation system is OSX or linux
if sys.platform in ["darwin", "win32"]:
	choice = input().lower()
elif sys.platform.startswith("linux"):
	choice = raw_input().lower()
else:
	print("unrecognized OS, check how to use raw input")
	choice = raw_input().lower()

import time
print("Start time:", time.ctime())
print("Finish time:", time.ctime())

######################
####     MAIN     ####
######################

from collections import OrderedDict
#creates an ordered dictionary
#but all dicts are ordered now!
seq_dict = OrderedDict()

#to avoid dictionary KeyErrors:
hightaxon = high_taxon_assignment_d.get(genus, "unassigned")

#modify df params:
def start():
	options = {
		"display": {
			"max_columns": None,
			"max_colwidth": 25,
			"expand_frame_repr": False,
			"max_rows": 14,
			"max_seq_items": 50,
			"precision": 4,
			"show_dimensions": True
		},
		"mode": {
			"chained_assignment": None
		}
	}
	for category, option in options.items():
		for op, value in option.items():
			pd.set_option(f"{category}.{op}", value)
if __name__ == "__main__":
	start()
	del start
#>> save as pandas_tricks.py, then call by export PYTHONSTARTUP="pandas_tricks.py" in Terminal?

#open a defined-size pandas dataframe with zeroes filled in:
df = pd.DataFrame(0.00, index=y_axis, columns=x_axis)
#fill iteratively
for x in some_list:
    for y in another_list:
        df.loc[x, y] = value

#make some fake data with test
import pandas.util.testing as tm
tm.N, tm.K = 15, 3 #rows, columns
tm.makeTimeDataFrame(freq="M").head()
[i for i in dir(tm) if i.startswith("make")]

#subset dataframe
sub = df[df["column"].isin(filterset)]
sub = df[df["value"] == filtervalue]

#join two existing dataframes on indexes
df = df1.join(df2, how="outer")
#join existing dataframes when the outer has multiple instances of the latter - autofill by station
bdf = bdf.join(sdf.set_index('Station'), on="Station#")
#alternatively, merge: 
df = df.merge(df1, how="left", left_on=["transcript_id"], right_on=["transcript_id"])

#write as tsv:
df.to_csv(path_or_buf='{}_maxdiffs_pd.out'.format(prefix), float_format='%0.3f', sep="\t")
df.to_pickle("./file.pkl.gz", compression="gzip")
df = pd.read_pickle("./file.pkl")
np.savetxt('{}_maxdiffs.out'.format(prefix), nparray, fmt='%0.3f', delimiter='\t')
#i.e. print tsv, with numbers truncated to three decimal places -> fmt/float_format

#to sort two lists in python:
ascendingtotals, ascendingtaxa = (list(x) for x in zip(*sorted(zip(unsortedList, unsortedTaxa)))) #parameters as reverse=True cannot be applied for *sorted

from Bio import SeqIO,AlignIO
#handling Fasta inputs
inFile = SeqIO.parse('input.txt', 'fasta')
from Bio.SeqUtils import GC
GC(rec.seq) #to calculate GC % of a sequence

with open('Split_fasta.txt', 'a') as result:
    for sequence in inFile:
        name = sequence.name
        seq = sequence.seq
        result.write('{}\t{}\n'.format(name,seq))
import gzip
with gzip.open('interproscan/' + f, "rt") as file:
	data = file.read() #from a gzipped file
with open("safe-seq.fasta") as infile, open("safe-seq.phy", "w") as outfile:
	alignmentfile = AlignIO.read(infile, "fasta") 
	#apparently also works as follows:
	# alignmentfile = AlignIO.read("some_file.fasta", "fasta")
	AlignIO.write(alignmentfile, outfile, "phylip-relaxed")
#for batch data:
files = os.listdir('.')
files = [f.replace(".fasta","") for f in files if f.endswith(".fasta")]
for file in files:
	with open("{}.fasta".format(file)) as infile, open("{}.phy".format(file), "w") as outfile:
		alignmentfile = AlignIO.read(infile, "fasta")
		AlignIO.write(alignmentfile, outfile, "phylip-relaxed")
alignmentfile.get_alignment_length() #gets alignment length
alignmentfile[:, x] #get column of this index, watch out this is 0-based

with open('a', 'w') as a, open('b', 'w') as b: # open multiple files at the same time


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#handling sequences
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()
protein = mySL.translate()



#to remove bad chars from sequence names
#and other useful things
badchars = (",:.;()'")
fullname = ''.join(c for c in fullname if c not in badchars)
#OR
prediction = re.sub('[():,]', '', prediction)


if sequence.startswith('>'):
	counter += 1


#to print shortened decimals
print('%.3f' % (somenumber))
print("{:.2f}".format(somenumber))


#avoid trailing newline in fastq
for index, item in enumerate(reads):
	if index != len(reads) - 1:
		name = item[0].split(' ')[0]
		out.write('{}/{}\n{}\n'.format(name, direction, '\n'.join(item[1:])))
	else:
		name = item[0].split(' ')[0]
		out.write('{}/{}\n{}\n+\n{}'.format(name, direction, item[1], item[3]))	


from tqdm import tqdm
for f in tqdm(files, desc="progress"):
	do whatever scripting

#or
with open(queryfile) as f:
	seqcount = f.read().count(">")
print("Sequence reading progress: {:.1f}%".format(100*(c / seqcount)))
 {:4.1f} #different notation - four characters altogether with the single decimal place

locals() #local variables?
dir(any_object) #print the structure of an object

############################
####  useful functions  ####
############################

#calculate centroid of a protein
from __future__ import division
import numpy as np

data = np.genfromtxt('file.pdb') #effective way how to import tables directly into an array
array = data[:, -3:] #if the data is here
np.mean(data[:,-3:], axis=0) #mean along the vertical axis, more or less means calculate centroid


def mean(lst):
	#mean without numpy
	if len(lst) > 0:
		return sum(lst) / len(lst)


def second(lst):
	#second-best in list
	return sorted(lst, reverse=True)[1]


def get_cmap(n, name='viridis'): #hsv for very divergent data?
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    colormap = plt.cm.get_cmap(name, n)
    rgbcolors = []
    for i in range(colormap.N):
        rgb = colormap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        rgbcolors.append(matplotlib.colors.rgb2hex(rgb))
    return rgbcolors
 
cmap = get_cmap(len(categories))
colors = {}
for i, X in enumerate(categories):
    colors[X] = cmap[i] #hopefully this is still recognized as a color
print(colors)
# also color = cmap.pop() works nicely



def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via raw_input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
		It must be "yes" (the default), "no" or None (meaning
		an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True,
			 "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '{}'" .format(default))

	while True:
		sys.stdout.write(question + prompt)
		if sys.platform in ["darwin", "win32"]:
			choice = input().lower()
		elif sys.platform.startswith("linux"):
			choice = raw_input().lower()
		else:
			print("unrecognized OS, check how to use raw input")
			choice = raw_input().lower()

		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")


#######################
####   FUNCTIONS   ####
#######################


#read only the first lines, without use of the counter
def read_ten(file_like_object):
    for line_number in range(10):
        x = file_like_object.readline()
        print(f"{line_number} = {x.strip()}")

###################
from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@gmail.com'

def force_taxid(accession):
	print("WARNING: missing taxid in input file taxified.x.out, requesting from NCBI server")
	prot = Entrez.efetch(db='protein', id=accession, rettype='gb', retmode='text')
	prot_record = SeqIO.read(prot, 'genbank')
	orgn = prot_record.annotations['organism']
	name2taxid = ncbi.get_name_translator([orgn])
	taxid = name2taxid[orgn][0]
	print("Organism retrieved:", orgn, taxid)

	return taxid


def delbadchars(string):
	""" remove unneeded characters """
	badchars = ("|+,:;()' ") #also []/@
	n = []
	for c in string:
		if c in badchars:
			c = "_"
		n.append(c)
	result = "".join(n)
	return result

def del_pattern(nodename):
	partpattern = r'p\d+'
	try:
		partID = re.search(partpattern, nodename).group()
		nodename = nodename.replace(partID, "")
		#this also works
		#for hit in re.findall(partpattern, nodename):
		#	nodename = nodename.replace(hit, "")
	except:
		pass
		#print(f"No pattern in {nodename}")
	return nodename

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def overlaps(moduleA, moduleB):
	startA = moduleA[0]
	endA = moduleA[1]
	startB = moduleB[0]
	endB = moduleB[1]
	if startA <= startB <= endA:
		if startA <= endB <= endA:
			overlap = "embed"
		else:
			overlap = "overlap"
	elif startA <= endB <= endA:
		overlap = "overlap"
	elif startB <= startA <= endB and startB <= endA <= endB:
		overlap = "embed"
	else:
		overlap = "none"
	return overlap

def overlaps(moduleA, moduleB):
	startA = moduleA[0]
	endA = moduleA[1]
	startB = moduleB[0]
	endB = moduleB[1]
	if startA <= startB <= endA or startA <= endB <= endA:
		overlap = "overlap"
	elif startB <= startA <= endB or startB <= endA <= endB:
		overlap = "overlap"
	else:
		overlap = "none"
	return overlap

def decor(string):
	def wrap():
		print("===============")
		print(string)
		print("===============")
	return wrap

def find_alignment_format(suffix):
	accepted = ["fasta", "phylip", "clustal", "emboss", "nexus", "stockholm"]
	if suffix in accepted:
		return suffix
	elif suffix in ["fasta", "fna", "faa", "fas", "fa"]:
		return "fasta"
	elif suffix in ["ali", "aln"]:
		return "fasta" #probably
	elif suffix in ["phylip", "phy"]:
		return "phylip-relaxed"
	elif suffix in ["nexus", "nex"]:
		return "nexus"
	else:
		quit("unrecognized MSA format")
