infile = open('Example_R1.fastq')
line = infile.read()
infile.close()

lines = line.split('\n')
nreads = len(lines) / 4

if len(lines) % 4 == 0:
	print "v poho"
else:
	print "neco je spatne"

reads = []

#fill the reads list
for i in range(nreads):
	#this will repeat the loop i times, which equals the number of sequences in fastq
	sequence = lines[:4]
	#vybere prvni ctyri polozky
	reads.append(sequence)
	lines = lines[4:]
	#odebere prvni ctyri polozky, ne 5
"""zkus tento prikaz:
l = range(16)
l[:4]
l[4:]
"""

out = open('Modified_R1.fastq', 'w')

#enumerate prida ke kazdemu item poradovy index, takhle muzu kontrolovat jestli uz nejsem na konci souboru
for index,item in enumerate(reads):
	if index != len(reads) - 1:
		newname = item[0].split(' ')[0]
		#tohle je zrychleny zapis vytvoreni promenne first = item.split(' ')[0] a newname = first[0]
		out.write(newname + "/1\n" + "\n".join(item[1:]) + "\n")
		#pozor na newline na konci, jeden radek pribyde a to by byl problem pro dalsi analyzy na tomto fastq vystupu
	else:
		newname = item[0].split(' ')[0] 
		out.write(newname + "/1\n" + "\n".join(item[1:]))
		#pozor na newline na konci, jeden radek pribyde a to by byl problem pro dalsi analyzy na tomto fastq vystupu


out.close()
