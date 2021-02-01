import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f, open(outfile, "w") as out:
	for line in f:
		out.write(line.replace("*", "X"))