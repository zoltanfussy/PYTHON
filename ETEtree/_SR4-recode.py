from Bio import AlignIO
#or maybe 

filename = "RESULT/SecY_v4_minusLBA_1-ali.fasta"

# groups = ["AGNPST","CHWY","DEKQR","FILMV"] #SR4
# recode = {}
# for i,aas in enumerate(groups):
# 	for aa in aas:
# 		recode[aa] = str(i+1)

recode = {'A': '1', 'G': '1', 'N': '1', 'P': '1', 'S': '1', 'T': '1', 
		'C': '2', 'H': '2', 'W': '2', 'Y': '2', 
		'D': '3', 'E': '3', 'K': '3', 'Q': '3', 'R': '3', 
		'F': '4', 'I': '4', 'L': '4', 'M': '4', 'V': '4',
		'-': '-', 'J': '-', 'X': '-'}


alignmentfile = AlignIO.read(filename, "fasta")
with open(filename.replace("-ali", "-recode"), "w") as result:
	for i,r in enumerate(alignmentfile):
		rname = r.id
		rseq = str(r.seq)
		for aa in recode:
			rseq = rseq.replace(aa, recode[aa])
		if i != len(alignmentfile) - 1:
			result.write(">{}\n{}\n".format(rname, rseq))
		else:
			result.write(">{}\n{}".format(rname, rseq))


