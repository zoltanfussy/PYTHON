infile = open('sequence.fasta')
line = infile.read()
infile.close()

#this is gonna create a very long string with all the file content on it
#we want to search for specific sequence in this string

seqs = line.split('>')
#the first item is empty, remove it with
seqs = seqs[1:]

infile = open('list1')
line = infile.read()
#we don`t need the line string anymore, it was saved as list "seqs"
infile.close()
names = line.split('\n')[:-1]
#in this input file, an extra newline is at the end of file, this will make the script ignore it

out = open('filtered_3.fasta','w')
#this is going to overwrite this file if existing, so careful

#declare an empty dictionary
seq_d = {}

#loop through sequences to fill the dictionary
for sequence in seqs:
	seqname =  sequence.split('\n')[0]
	seq_d[seqname] = sequence

for name in names:
 	print seq_d[name]

out.close()