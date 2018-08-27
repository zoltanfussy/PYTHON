#create function to read in fastq files based on rename_fastq-range.py
#do it so that you can input file name and new appendix to the name: "/1" or "/2"
#not too handy considering that it loads the entire file into memory...
def fastq_fill_list(lines,nreads):
	print "fastq_fill_list called by fastq_modify.\n..."
	reads = []
	for i in range(nreads):
		sequence = lines[:4]
		reads.append(sequence)
		lines = lines[4:]
	print "Output returned by fastq_fill_list."
	return reads

def fastq_modify(infile,tag):
	print "This script modifies input FASTQ by adding a tag at the end of nameline.\n"
	infilename = infile
	print "Reading file " + infilename + ", tagging by " + tag + ".\n..."
	infile = open(infile)
	line = infile.read()
	infile.close()

	lines = line.split('\n')
	nreads = len(lines) / 4
	reads = fastq_fill_list(lines,nreads)

	outfilename = "Modified_" + infilename.replace('.fastq','')[0] + ".fastq"
	out = open(outfilename, 'w')

	for index, item in enumerate(reads):
			if index != len(reads) - 1:
				newname = item[0].split(' ')[0]
				out.write(newname + tag + "\n" + "\n".join(item[1:]) + "\n")
			else:
				newname = item[0].split(' ')[0] 
				out.write(newname + tag + "\n" + "\n".join(item[1:]))

	print "Modified sequences written to " + outfilename + ".\nThank you for using fastq_modify. God bless you in abundance.\n\n\n---"
	out.close()



fastq_modify('Example_R1.fastq','/1')
fastq_modify('Example_R2.fastq','/2')