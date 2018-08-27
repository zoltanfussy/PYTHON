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

out = open('filtered.fasta','w')
#this is going to overwrite this file if existing, so careful

# loop through sequences 
for sequence in seqs:
	seqname =  sequence.split('\n')[0]
	#loop through query names
	for item in names:
		if item == seqname:
			out.write ('>' + sequence)
			#print '>' + sequence
		else:
			pass
# tohle by mohlo taky fungovat :D
#	if seqname in names:
#		aminoacids = ''.join(sequence[1:])
#		print ('>') + seqname + "\n" + aminoacids
out.close()
#but this is very slow on thousands of sequences - use dictionaries for that