#this script renames PASTA alignments back to original, but takes only 30 character substring of the name! mainly useful for my house-made BLAST database
import argparse

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--prefix', help='PASTA prefix', required=True)

args = parser.parse_args()

prefix = args.prefix

print("PASTA prefix defined:%s" %(prefix))
intranslation = open(prefix + '_temp_name_translation.txt')
line = intranslation.read()
intranslation.close()

codes = line.split('\n')
nreads = int(len(codes) // 3)

translations = {}

#fill the reads list
for i in range(nreads):
	#this will repeat the loop i times, which equals the number of sequences in the dataset
	#only first 30 chars are taken so the script is compatible with the prediction tree renamer
	translations[codes[0]] = codes[1][:30]
	#removes first three lines
	codes = codes[3:]

print("Reading alignment to rename...")
inalign = open(prefix + '_temp_iteration_2_seq_alignment.txt')
line = inalign.read()
inalign.close()

lines = line.split('\n')

if len(lines) % 2 == 0:
	pass
else:
	lines = lines[:-1]


outfile = open(prefix + '.alignment.fasta', 'w')

for sequence in lines:
	if sequence[0] == '>':
		safeheader = sequence[1:]
		header = translations[safeheader]
		outfile.write('>' + header + '\n')
	else:
		outfile.write(sequence + '\n')

outfile.close()
print("Alignment written to {}.alignment.fasta. \n Reading tree to rename...".format(prefix))

strom = open(prefix + '_temp_iteration_2_tree.tre')
strom_line = strom.readline()

for key in translations:
    strom_line = strom_line.replace(key, translations[key])


with open(prefix + '.some.tre', 'w') as result:
    result.write(strom_line)

print("Tree written to {}.some.tree. Terminating.".format(prefix))