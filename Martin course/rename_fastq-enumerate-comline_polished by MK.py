import argparse

parser = argparse.ArgumentParser(description='Modify fastq file')
parser.add_argument('-i', '--infile', help='input filename', required=True)
parser.add_argument('-p', '--pair', help='specifies pair 1 or 2', required=True)
parser.add_argument('-o', '--outfile', help='output filename', required=True)
args = parser.parse_args()
infile = open(args.infile)
line = infile.read()
infile.close()

lines = line.split('\n')
nlines = len(lines)
number_reads = nlines/4


reads = []

for i in range(number_reads):
	sequence = lines[:4]
	reads.append(sequence)
	lines = lines[4:]
	
print reads
print len(reads)
direction = str(args.pair)
out = open(args.outfile, 'w')
 
for index, item in enumerate(reads):
	if index != len(reads) - 1:
		name = item[0].split(' ')[0]
		out.write('%s/%s\n%s\n' % (name, direction, '\n'.join(item[1:])))
	else:
		name = item[0].split(' ')[0]
		out.write('%s/%s\n%s\n+\n%s' % (name, direction, item[1], item[3]))		
out.close()

