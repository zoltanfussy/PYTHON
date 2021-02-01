import argparse, gzip, time

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print("starting {}".format(current_time))

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Diamond outfile to be taxified', required=True)
args = parser.parse_args()
infile = args.infile

with gzip.open(infile, mode='rb') as fastqfile, open(infile.replace("fastq.gz", "fasta"), "w") as subset:
	print("reading fastq...")
	c = 0
	for l in fastqfile:
		c += 1
		if c % 10000000 == 0:
			print("{}M".format(c/1000000))
		l = l.strip().decode('utf8')
		try:
			if c % 4 == 1:
				l = l.split()[0]
				subset.write("{}\n".format(l.replace("@", ">")))
			elif c % 4 == 2:
				subset.write("{}\n".format(l))
		except IndexError:
			print(l, "not enough columns")