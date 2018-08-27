infile1 = open('Example_R1.fastq')
line1 = infile1.read()
infile1.close()

out1 = open('processed_R1.fastq','w')

seqs1 = line1.split('\n')
for lines in seqs1:
	if lines.startswith(@SRR):
		newhandle = '.'.join(lines.split()[0].split('.')[0-1])
		out1.write(newhandle + "/1\n")
	else:
		out1.write(lines + "\n")

#to same pro druhy soubor 
infile2 = open('Example_R2.fastq')
line2 = infile2.read()
infile2.close()

out2 = open('processed_R2.fastq','w')

seqs2 = line2.split('\n')
for lines in seqs2:
	if lines.startswith("@MISEQ"):
		newhandle = lines.split(' ')[0]
		out2.write(newhandle + "/2\n")
		#pozor, kdybys chtel vytisknout backslash, tak musis napsat "\\2\n"
	else:
		out2.write(lines + "\n")

out1.close()
out2.close()

#jiny pristup je spojit jednotlive polozky fastq do jedne item v seznamu pomoci split("@MISEQ"), 
#dale for loop s nasledujicimi kroky
#1) extrahovat jmeno sekvence bez konce - to co je za mezerou
#2) spojit zbytek item pomoci '\n'.join(sekvence.split('\n')[1])
#3) vytisknout jmeno + "\n" + zbytek + "\n"