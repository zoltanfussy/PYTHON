# Performs fasta header renaming, and calls mafft and trimal to align and trim datasets 
# automatically on all fasta files in the current folder
# Alternatively, set folder with -d

import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-d', '--directory', help='Working directory', default='.')
parser.add_argument('-s', '--skip', help='File with sequence names to skip (names will be processed with split("_")', default='')
args = parser.parse_args()
directory = args.directory

#### Functions ####
###################

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
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

###############################
#### Open and parse inputs ####
###############################

files = os.listdir(directory)
files = [f for f in files if f.endswith(".fasta")]

filecount = len(files)


#eventually, remove bad datasets and bad chars from names
if args.skip == '':
	skip = ['cilMESOD', 'redCYANm']
else:
	skip = []
	with open(args.skip) as inskip:
		for line in inskip:
			skip.append(line.strip().split('_')[0])
#badchars = (",:.;()'")
#fullname = ''.join(c for c in fullname if c not in badchars) #place this somewhere in the script
print(", ".join(skip))
if query_yes_no("Is skiplist okay?") is True:
	print("continuing...")
	pass
else:
	quit("script stopped!")

fnames = set()
for f in files:
	leaves = []
	fname = f.split(".fasta")[0]
	fnames.add(fname)
	with open(fname + "-safe.fst", "w") as outsafefile, open(fname + "-trans.tsv", "w") as outtransfile:
		for sequence in SeqIO.parse(f, "fasta"):
			safeleaf = str(sequence.name).split("_")[0]
			leaves.append(safeleaf)
			if leaves.count(safeleaf) > 1:
				print("Duplicate found: " + safeleaf)
				count = leaves.count(safeleaf)
				safeleaf = safeleaf + str(count)
			if safeleaf not in skip:
				outsafefile.write(">{}\n{}\n".format(safeleaf, sequence.seq))
				outtransfile.write("{}\t{}\n".format(safeleaf, sequence.description))
	print("File {} processed".format(f))

print("Writing sequences done on {} files, now to alignment and trimming...".format(filecount))

zerosizefiles = []
for f in fnames:
	fname = f + "-safe.fst"
	if os.stat(fname).st_size > 0:
		os.system("mafft --maxiterate 1000 --localpair --thread 3 {0}-safe.fst > {0}.linsi".format(f))
		os.system("trimal -in {0}.linsi -out {0}.MaffTrimal.fst -fasta -automated1".format(f))
	else:
		zerosizefiles.append(fname)

for fname in zerosizefiles:
	if os.stat(fname).st_size > 0:
		os.system("mafft --maxiterate 1000 --localpair --thread 3 {0}-safe.fst > {0}.linsi".format(fname))
		os.system("trimal -in {0}.linsi -out {0}.MaffTrimal.fst -fasta -automated1".format(fname))
	else:
		print("Could not process {}. Zero bytes in file".format(fname))
